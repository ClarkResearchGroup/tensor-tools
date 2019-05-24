//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef AUTOMPO
#define AUTOMPO

#include <algorithm>
#include <map>
#include <type_traits>
#include "AutoMPO.h"

using std::find;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::array;
using std::pair;
using std::make_pair;
using std::move;
using std::min;
using std::max;
using std::map;
using std::set;

using Cplx = std::complex<double>;
using Real = double;
void Error(string s){
  perr<<s<<endl;
  assert(1==2);
}

bool
isZero(Cplx const& z, Real thresh = 1E-13) { return std::abs(z) < thresh; }

bool
less(Real x, Real y, Real eps = 1E-12)
    {
    Real ax = std::fabs(x);
    Real ay = std::fabs(y);
    Real scale = (ax < ay ? ay : ax);
    return (y-x) >= scale*eps;
    }

bool
equal(Cplx x, Cplx y, Real eps = 1E-12)
    {
    Real ax = std::abs(x);
    Real ay = std::abs(y);
    Real scale = (ax < ay ? ay : ax);
    return std::abs(x-y) <= scale*eps;
    }

bool
less(Cplx const& z1, Cplx const& z2, Real eps = 1E-12) 
    { 
    if(not equal(z1.real(),z2.real(),eps)) 
        {
        return less(z1.real(),z2.real(),eps);
        }
    return less(z1.imag(),z2.imag());
    }

bool
isReal(const Cplx& z) { return z.imag() == 0; }

bool
isApproxReal(Cplx const& z, Real epsilon = 1E-12) { return std::fabs(z.imag()) < epsilon; }

SiteTerm::
SiteTerm() : i(-1) { }

SiteTerm::
SiteTerm(string const& op_,
         int i_)
    :
    op(op_),
    i(i_)
    { }

bool
isFermionic(SiteTerm const& st)
    {
#ifdef DEBUG
    for(char c : st.op)
    if(c == '*')
        {
          std::cerr<<st.op<<std::endl;
        Error("SiteTerm contains a '*' but isFermionic does not handle this case");
        }
#endif
    if(!st.op.empty() && st.op.front() == 'C') return true;
    return false;
    }

SiteTermProd 
mult(SiteTermProd first, 
     SiteTermProd const& second)
    {
    first.insert(first.end(), second.begin(), second.end());
    return first;
    }    

bool
isFermionic(SiteTermProd const& sprod)
    {
    bool isf = false;
    for(auto& st : sprod)
        {
        //Flip isf in a Z2 fashion for every fermionic operator
        if(isFermionic(st)) isf = !isf;
        }
    return isf;
    }

string
fermionicTerm(const string& op)
    {
    static array<pair<string,string>,6>
           rewrites =
           {{
           make_pair("Cdagup","Adagup"),
           make_pair("Cup","Aup"),
           make_pair("Cdagdn","Adagdn*Fup"),
           make_pair("Cdn","Adn*Fup"),
           make_pair("C","A"),
           make_pair("Cdag","Adag")
           }};
    for(auto& p : rewrites)
        {
        if(p.first == op) return p.second;
        }
    return op;
    }

void 
rewriteFermionic(SiteTermProd & prod, 
                 bool isleftFermionic)
    {
    if(prod.empty()) Error("Empty product in rewriteFermionic is not expected.");    
    
    int i = prod.front().i;
    for(auto& st : prod)
        if(st.i != i)
            {
            Error("Multi-site product in rewriteFermionic is not expected.");    
            }

    // Rewrite a fermionic single site product using the Jordan-Wigner string            
    bool isSiteFermionic = isFermionic(prod);
    if(isSiteFermionic)
        {
        for(auto& st : prod) if(isFermionic(st)) st.op = fermionicTerm(st.op);
        }
    
    // Add a FermiPhase operator at the end if the product of operators
    // to the left (including this site) is fermionic
    if((isleftFermionic && !isSiteFermionic) || (!isleftFermionic && isSiteFermionic))
        {
        prod.emplace_back("F", i);         
        }
    }


void HTerm::
add(string const& op,
    int i,
    Real x)
    {
    //The following ensures operators remain
    //in site order within the vector "ops"
    auto it = ops.begin();
    while(it != ops.end() && it->i <= i) ++it;

    auto t = SiteTerm(op,i);

    // If the operator is fermionic and being inserted in between existing operators 
    // need to check if an extra minus is required
    if(it != ops.end() && isFermionic(t))
        { 
        auto rightOps = SiteTermProd(it,ops.end());
        if(isFermionic(rightOps)) coef *= -1;
        }   

    coef *= x;
    ops.insert(it,t);
    }

HTerm& HTerm::
operator*=(Real x)
    {
    if(Nops() == 0) Error("No operators in HTerm");
    coef *= x;
    return *this;
    }

HTerm& HTerm::
operator*=(Cplx x)
    {
    if(Nops() == 0) Error("No operators in HTerm");
    coef *= x;
    return *this;
    }

bool HTerm::
operator==(HTerm const& o) const
    {
    if(not equal(coef,o.coef,1E-12)) return false;
    if(Nops() != o.Nops()) return false;

    for(size_t n = 0; n <= ops.size(); ++n)
    if(ops[n] != o.ops.at(n)) 
        {
        return false;
        }

    return true;
    }

bool HTerm::
operator<(HTerm const& o) const 
    { 
    if(not equal(coef,o.coef,1E-12)) return less(coef,o.coef,1E-12);
    return (ops < o.ops);
    }

bool LessNoCoef::
operator()(HTerm const& t1, HTerm const& t2) const
    { 
    if(t1.ops.size() != t2.ops.size()) return t1.ops.size() < t2.ops.size();
            
    for(size_t j = 0ul; j < t1.ops.size(); ++j)
        {
        if(t1.ops[j] != t2.ops[j]) return t1.ops[j] < t2.ops[j];
        }
    return false;
    }

AutoMPO::Accumulator::
Accumulator(AutoMPO* pa_, 
            Real x_)
    :
    pa(pa_),
    state(New),
    coef(x_)
    {}

AutoMPO::Accumulator::
Accumulator(AutoMPO* pa_, 
            Cplx x_)
    :
    pa(pa_),
    state(New),
    coef(x_)
    {}

AutoMPO::Accumulator::
Accumulator(AutoMPO* pa_)
    : 
    Accumulator(pa_,1)
    {}


AutoMPO::Accumulator::
Accumulator(AutoMPO* pa_, 
            const char* op_)
    :
    pa(pa_),
    state(Op),
    coef(1),
    op(op_)
    {}

AutoMPO::Accumulator::
Accumulator(AutoMPO* pa_, 
            const string& op_)
    :
    pa(pa_),
    state(Op),
    coef(1),
    op(op_)
    {}


AutoMPO::Accumulator::
~Accumulator()
    {
    if(state==Op) Error("Invalid input to AutoMPO (missing site number?)");
    term *= coef;
    pa->add(term);
    }
    

AutoMPO::Accumulator& AutoMPO::Accumulator::
operator,(Real x)
    {
    coef *= x;
    return *this;
    }

AutoMPO::Accumulator& AutoMPO::Accumulator::
operator,(Cplx x)
    {
    coef *= x;
    return *this;
    }

AutoMPO::Accumulator& AutoMPO::Accumulator::
operator,(int i)
    {
    if(state==Op)
        {
        term.add(op,i);
        state = New;
        op = "";
        }
    else
        {
        coef *= Real(i);
        }
    return *this;
    }

AutoMPO::Accumulator& AutoMPO::Accumulator::
operator,(const char* op_)
    {
    if(state == New)
        {
        op = op_;
        state = Op;
        }
    else
        {
        Error("Invalid input to AutoMPO (two strings in a row?)");
        }
    return *this;
    }

AutoMPO::Accumulator& AutoMPO::Accumulator::
operator,(const string& op_)
    {
    if(state == New)
        {
        op = op_;
        state = Op;
        }
    else
        {
        Error("Invalid input to AutoMPO (two strings in a row?)");
        }
    return *this;
    }

void AutoMPO::
add(HTerm & t)
    {
    if(abs(t.coef) == 0.0) return;

    //go through and add identities
    int l = t.ops.front().i;
    for(int idx=0;idx<t.ops.size();idx++){
      auto currSpot = t.ops[idx].i;
      if(l<currSpot){
        while(l!=currSpot)
          t.add("Id",l++);
      }
      else
      l++;
    }
    auto it = terms_.find(t);
    if(it == terms_.end())
        {
        terms_.insert(move(t));
        }
    else //found duplicate
        {
        auto nt = t;
        nt.coef += it->coef;
        terms_.erase(it);
        terms_.insert(move(nt));
        }
    }

/*
MPO convention:
===============
For each link of the MPO, define a set of bases 
that describe the terms of the Hamiltonian
corresponding to the left "half" of the MPO.
The terms include "IL", which means the product
of identities to the left, and "HL", the sum of
all terms entirely contained on the left.

Usually these two special terms occupy positions 1 and two,
respectively.

The rest of the bases are each site term on the left that
is connected to something on the right.

So for neighbor and next neighbor, operator pair A B, 
coefs t1 and t2, on site n, the MPO matrix is:
n-1             n
      1111   HL  11A1  111A  <== bases
1111   1     0     0    A         
HL     0     1     0    0   
11A1   0    t2 B   0    0   
111A   0    t1 B   1    0   

For neighbor and next neighbor, operator pair A B and B A, t1 and t2
site n:
n-1             n
      1111  HL    11A1 11B1  111A  111B 
1111   1     0     0     0     A     B  
HL     0     1     0     0     0     0  
11A1   0    t2 B   0     0     0     0  
11B1   0    t2 A   0     0     0     0  
111A   0    t1 B   1     0     0     0  
111B   0    t1 A   0     1     0     0  

F == fermiPhase, i.e. F = (-1)^(# of fermions of either type of spin)
Then we make c and cdagger both have F's going off to the left.

Fermion operator rewriting convention:

//
//Spinless fermions
//

Cdag_i C_j  = (F_1 F_2 F_3 ... F_{i-1})^2 (Adag_i F_i) F_{i+1} ... A_j
            = Adag_i F_{i+1} ... A_j

C_i Cdag_j = (A_i F_i) F_{i+1} ... Adag_j

//
//Fermions with spin
//

Cdagup_i Cup_j  = (F_1 F_2 F_3 ... )^2 (Adagup_i F_i) F_{i+1} ... Aup_j
                = (Adagup_i F_i) F_{i+1} ... Aup_j //cancel squared F operators

Cup_i Cdagup_j = (Aup_i F_i) F_{i+1} ... Adagup_j

Cdagdn_i Cdn_j  = (Adagdn_i F_i) F_{i+1} ... Fup_j Adn_j 
                = - Adagdn_i F_{i+1} ... Fup_j Adn_j     //use Adagdn_i * F_i = -Adagdn_i
                = Adagdn_i F_{i+1} ... Fup_j Fdn_j Adn_j //use Adn_j = -Fdn_j*Adn_j
                = Adagdn_i F_{i+1} ... (F_j Adn_j)       //combine Fup_j*Fdn_j = F_j (definition)

Cdn_i Cdagdn_j = (Adn_i F_i) F_{i+1} ... Fup_j Adagdn_j
               = - Adn_i F_{i+1} ... Fup_j Adagdn_j      //use Adn_i*F_i = -Adn_i
               = Adn_i F_{i+1} ... Fup_j Fdn_j Adagdn_j  //use Adagdn_j = -Fdn_j*Adagdn_j
               = Adn_i F_{i+1} ... (F_j Adagdn_j)        //combined Fup_j*Fdn_j = F_j (definition)


*/


//TODO:
// o Add support for > 2 site operators
// o Add support for long-range (exponentially-decaying type) operator strings
// o Add support for fermionic operator strings


void
plusAppend(string & s, string const& a)
    {
    if(s.size() == 0 || s == "0") s = a;
    else 
        {
        s += "+";
        s += a;
        }
    }

//#define SHOW_AUTOMPO


string
startTerm(const string& op)
    {
    static array<pair<string,string>,6>
           rewrites =
           {{
           make_pair("Cdagup","Adagup*F"),
           make_pair("Cup","Aup*F"),
           make_pair("Cdagdn","Adagdn"),
           make_pair("Cdn","Adn"),
           make_pair("C","A*F"), //A*F is -A, so essentially a trick for putting in a -1
           make_pair("Cdag","Adag")
           }};
    for(auto& p : rewrites)
        {
        if(p.first == op) return p.second;
        }
    return op;
    }

string
endTerm(const string& op)
    {
    static array<pair<string,string>,6>
           rewrites =
           {{
           make_pair("Cup","Aup"),
           make_pair("Cdagup","Adagup"),
           make_pair("Cdn","F*Adn"),
           make_pair("Cdagdn","F*Adagdn"),
           make_pair("C","A"),
           make_pair("Cdag","Adag")
           }};
    for(auto& p : rewrites)
        {
        if(p.first == op) return p.second;
        }
    return op;
    }


//

void 
decomposeTerm(int n, 
              SiteTermProd const& ops, 
              SiteTermProd & left, 
              SiteTermProd & onsite, 
              SiteTermProd & right)
    {
    auto isOnSiteOrOnTheRight = [&n](const SiteTerm &t) {return t.i >= n;};
    auto startOfOnSite = find_if(ops.begin(), ops.end(), isOnSiteOrOnTheRight);
    
    auto isOnTheRight = [&n](const SiteTerm &t) {return t.i > n;};
    auto startOfRightPart = find_if(startOfOnSite, ops.end(), isOnTheRight);

    left = SiteTermProd(ops.begin(), startOfOnSite);
    onsite = SiteTermProd(startOfOnSite, startOfRightPart);
    right = SiteTermProd(startOfRightPart, ops.end());
    }  




/*template<typename T, typename V>
auto
forceType(V x) -> typename std::enable_if<std::is_same<T,V>::value,T>::type
    {
    return x;
    }
template<typename T,
         class = std::enable_if<std::is_same<T,double>::value>,T>::type
T
forceType(Cplx z) { return z.real(); }*/

Real
conj(Real x) { return x; }


std::ostream& 
operator<<(std::ostream& s, SiteTerm const& t)
    {
    s << t.op << "(" << t.i << ")";
    return s;
    }


std::ostream& 
operator<<(std::ostream& s, HTerm const& t)
    {
    const char* pfix = "";
    if(abs(t.coef-1.0) > 1E-12) 
        { //HACK
        s << t.coef.real(); //(isReal(t.coef) ? t.coef.real() : t.coef);
        }
    for(auto& st : t.ops) 
        {
        s << pfix<< st.op<< "("<<st.i<<")"; //printf("%s%s(%d)",pfix,st.op,st.i);
        pfix = " ";
        }
    return s;
    }

std::ostream& 
operator<<(std::ostream& s, const AutoMPO& a)
    {
    s << "AutoMPO:\n";
    for(const auto& t : a.terms()) s << t << "\n";
    return s;
    }

#endif
