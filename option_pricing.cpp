// Option Pricing Lab: a portfolio refactor of Baruch C++ coursework.
// See README.md for provenance, model assumptions, and validation.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace quant {
enum class Type { Call, Put };
struct Option {
    double spot, strike, rate, carry, volatility, maturity;
    Type type;
    void validate() const {
        for (double x : {spot, strike, rate, carry, volatility, maturity})
            if (!std::isfinite(x)) throw std::invalid_argument("Parameters must be finite");
        if (spot <= 0 || strike <= 0 || volatility < 0 || maturity < 0)
            throw std::invalid_argument("Spot/strike must be positive; volatility/maturity nonnegative");
    }
    double payoff(double s) const {
        return std::max((type == Type::Call ? 1.0 : -1.0) * (s - strike), 0.0);
    }
};
double normal_cdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }
struct Analytics { double price, delta, gamma; };
Analytics black_scholes(const Option& o) {
    o.validate();
    const double growth = std::exp((o.carry-o.rate)*o.maturity);
    const double discount = std::exp(-o.rate*o.maturity);
    const double forward = o.spot * std::exp(o.carry*o.maturity);
    const double sign = o.type == Type::Call ? 1.0 : -1.0;
    if (o.maturity == 0 || o.volatility == 0) {
        // At the payoff kink, delta/gamma are undefined; report NaN.
        if (forward == o.strike) return {0, NAN, NAN};
        return {discount*o.payoff(forward), sign*(sign*(forward-o.strike)>0 ? growth : 0), 0};
    }
    const double v = o.volatility*std::sqrt(o.maturity);
    const double d1 = (std::log(o.spot/o.strike)+(o.carry+0.5*o.volatility*o.volatility)*o.maturity)/v;
    const double d2 = d1-v;
    const double density = std::exp(-0.5*d1*d1)/std::sqrt(2*std::acos(-1.0));
    return {sign*(o.spot*growth*normal_cdf(sign*d1)-o.strike*discount*normal_cdf(sign*d2)),
            sign*growth*normal_cdf(sign*d1), growth*density/(o.spot*v)};
}
struct Estimate { double price, standard_deviation, standard_error; std::size_t pairs; };
Estimate monte_carlo(const Option& o, std::size_t pairs, std::uint64_t seed) {
    o.validate();
    if (pairs < 2) throw std::invalid_argument("At least two independent pairs required");
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal;
    const double drift=(o.carry-0.5*o.volatility*o.volatility)*o.maturity;
    const double diffusion=o.volatility*std::sqrt(o.maturity);
    const double discount=std::exp(-o.rate*o.maturity);
    double mean=0, m2=0;
    for (std::size_t i=0; i<pairs; ++i) {
        const double z=normal(rng);
        // Treat each antithetic pair as ONE independent observation for the SE.
        const double value=0.5*discount*(o.payoff(o.spot*std::exp(drift+diffusion*z))+
                                                   o.payoff(o.spot*std::exp(drift-diffusion*z)));
        const double difference=value-mean;
        mean+=difference/static_cast<double>(i+1);
        m2+=difference*(value-mean); // Welford: constant memory, stable variance.
    }
    const double sd=std::sqrt(m2/static_cast<double>(pairs-1));
    return {mean, sd, sd/std::sqrt(static_cast<double>(pairs)), pairs};
}
// Explicit time stepping with upwind drift: nonnegative stencil weights.
// Automatic time-step selection enforces monotonicity for either carry sign.
double finite_difference(const Option& o, std::size_t cells=300) {
    o.validate();
    if (cells < 10 || cells > 2000) throw std::invalid_argument("Grid cells must be 10..2000");
    if (o.maturity==0 || o.volatility==0) return black_scholes(o).price;
    const double upper=4*std::max(o.spot,o.strike), ds=upper/static_cast<double>(cells);
    const double m=static_cast<double>(cells);
    const double bound=o.volatility*o.volatility*m*m+std::abs(o.carry)*m+std::max(o.rate,0.0);
    const double required=std::ceil(o.maturity*bound/0.9);
    if (!std::isfinite(required) || required>2000000) throw std::invalid_argument("PDE grid exceeds work budget");
    const auto steps=static_cast<std::size_t>(std::max(1.0,required));
    const double dt=o.maturity/static_cast<double>(steps);
    std::vector<double> old(cells+1), next(cells+1);
    for (std::size_t i=0;i<=cells;++i) old[i]=o.payoff(static_cast<double>(i)*ds);
    for (std::size_t n=1;n<=steps;++n) {
        const double tau=static_cast<double>(n)*dt;
        next[0]=o.type==Type::Put ? o.strike*std::exp(-o.rate*tau) : 0;
        next[cells]=o.type==Type::Call ? std::max(upper*std::exp((o.carry-o.rate)*tau)-o.strike*std::exp(-o.rate*tau),0.0) : 0;
        for (std::size_t i=1;i<cells;++i) {
            const double j=static_cast<double>(i), diffusion=0.5*o.volatility*o.volatility*j*j;
            const double lower=diffusion+std::max(-o.carry*j,0.0);
            const double higher=diffusion+std::max(o.carry*j,0.0);
            next[i]=dt*lower*old[i-1]+(1-dt*(lower+higher+o.rate))*old[i]+dt*higher*old[i+1];
        }
        old.swap(next);
    }
    const double position=o.spot/ds;
    const auto i=static_cast<std::size_t>(position);
    return old[i]+(position-static_cast<double>(i))*(old[i+1]-old[i]);
}
double perpetual_american(const Option& o) {
    o.validate();
    if (o.rate<=0 || o.volatility<=0 || o.carry>=o.rate)
        throw std::invalid_argument("Perpetual engine requires r>0, sigma>0, b<r");
    const double v2=o.volatility*o.volatility, a=0.5-o.carry/v2;
    const double root=std::sqrt(a*a+2*o.rate/v2);
    const double y=o.type==Type::Call ? a+root : a-root;
    const double boundary=o.strike*y/(y-1);
    const bool exercise=o.type==Type::Call ? o.spot>=boundary : o.spot<=boundary;
    return exercise ? o.payoff(o.spot) : o.payoff(boundary)*std::pow(o.spot/boundary,y);
}
void check(bool pass, const std::string& name) {
    if (!pass) throw std::runtime_error("Test failed: "+name);
}
void tests() {
    Option o{100,100,0.05,0.05,0.2,1,Type::Call};
    check(std::abs(black_scholes(o).price-10.4505835721856)<1e-10,"reference call");
    for (Option x : std::vector<Option>{{60,65,.08,.08,.3,.25,Type::Call},o,
                          {102,122,.045,0,.43,1.65,Type::Call},
                          {100,100,.03,-.02,.2,1,Type::Call}}) {
        const auto c=black_scholes(x); x.type=Type::Put; const auto p=black_scholes(x);
        check(std::abs(c.price-p.price-x.spot*std::exp((x.carry-x.rate)*x.maturity)+x.strike*std::exp(-x.rate*x.maturity))<1e-10,"generalized parity");
        const double h=.01; auto up=x, down=x; up.spot+=h; down.spot-=h;
        check(std::abs((black_scholes(up).price-black_scholes(down).price)/(2*h)-p.delta)<1e-6,"delta finite difference");
        check(std::abs((black_scholes(up).price-2*p.price+black_scholes(down).price)/(h*h)-p.gamma)<1e-6,"gamma finite difference");
        for (Type type : {Type::Call,Type::Put}) {
            x.type=type; const double exact=black_scholes(x).price;
            check(std::abs(finite_difference(x,400)-exact)<.10,"PDE benchmark");
            const auto mc=monte_carlo(x,100000,42);
            check(std::abs(mc.price-exact)<6*mc.standard_error+.001,"MC statistical benchmark");
        }
    }
    auto edge=o; edge.maturity=0;
    check(black_scholes(edge).price==0 && std::isnan(black_scholes(edge).delta),"expiry kink");
    edge=o; edge.volatility=0;
    check(std::abs(black_scholes(edge).price-(100-100*std::exp(-.05)))<1e-10,"zero volatility");
    bool rejected=false; edge.spot=-1;
    try { (void)black_scholes(edge); } catch(const std::invalid_argument&) { rejected=true; }
    check(rejected,"invalid spot");
    Option american{110,100,.1,.02,.1,1,Type::Call};
    check(std::abs(perpetual_american(american)-18.5034998830479)<1e-8,"perpetual reference");
    american.spot=1000; check(perpetual_american(american)==900,"call exercise region");
    american.type=Type::Put; american.spot=1; check(perpetual_american(american)==99,"put exercise region");
    const double coarse=std::abs(finite_difference(o,100)-black_scholes(o).price);
    const double fine=std::abs(finite_difference(o,400)-black_scholes(o).price);
    check(fine<coarse,"PDE refinement");
    std::cout<<"All validation checks passed.\n";
}
void demo() {
    std::cout<<std::fixed<<std::setprecision(6);
    std::cout<<"type,spot,strike,rate,carry,volatility,maturity,analytic,delta,gamma,mc,mc_se,ci95_low,ci95_high,pde\n";
    for (double spot : {80.0,100.0,120.0}) for (Type type : {Type::Call,Type::Put}) {
        Option o{spot,100,.05,.03,.2,1,type};
        const auto a=black_scholes(o); const auto mc=monte_carlo(o,100000,42);
        std::cout<<(type==Type::Call ? "call" : "put")<<','<<spot<<",100,0.05,0.03,0.2,1,"<<a.price<<','<<a.delta<<','<<a.gamma<<','<<mc.price<<','<<mc.standard_error<<','<<mc.price-1.96*mc.standard_error<<','<<mc.price+1.96*mc.standard_error<<','<<finite_difference(o)<<'\n';
    }
}
} // namespace quant
int main(int argc,char** argv) {
    try {
        const std::string mode=argc==1 ? "--demo" : argv[1];
        if (argc>2) throw std::invalid_argument("Expected at most one argument");
        if (mode=="--test") quant::tests();
        else if (mode=="--demo") quant::demo();
        else if (mode=="--help") std::cout<<"Usage: option_pricing [--demo|--test|--help]\nDemo emits CSV. Parameters are defined in quant::demo().\n";
        else throw std::invalid_argument("Unknown option; use --help");
        return 0;
    } catch(const std::exception& e) { std::cerr<<"Error: "<<e.what()<<'\n'; return 1; }
}
