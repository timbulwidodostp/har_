program define kernel_F_all, rclass

syntax anything [if][in],kernel(string)[,noConstant Level(cilevel)]
marksample touse


qui{

	local opnames B P Q 
	local ops B P Q
	
	
	if "`kernel'"==""{
    display as err "Error:kernel must be input"
	error 198
	}
	else {
	local nop : list posof "`kernel'" in opnames
	if !`nop'{
	display as err "Error:kernel must be chosen from `opnames'"
	error 198
	}
	local op :word `nop' of `ops'
	}

	
	_iv_parse `0'

	local lhs `s(lhs)'
	local endog `s(endog)'
	local exog `s(exog)'
	local inst `s(inst)'

    
	markout `touse' `lhs' `exog' `inst' `endog'


   local endo_ct : word count `endog'
   local exog_ct : word count `exog'
   local inst_ct : word count `inst'

    if `endo_ct' > `inst_ct' {
		di as err "equation not identified; must have at " "least as many instruments not in"
		di as err "the regression as there are "  "instrumented variables"
		exit 481
	}
	if (`endo_ct' + `exog_ct' == 0) & "`noconstant'"!="" {
		di as err "no RHS variables"
		exit 498
	}
	

	if `"`constant'"'==`""' {
	tempvar CONS
	gen byte `CONS' = 1
	local carg `""'
	}
	else {
	local CONS `""'
	local carg `"nocons"'
    }
	
	
	local depva `lhs'
	local indepv `endog' `exog' `CONS'
	local instv  `inst' `exog' `CONS'
	
	
local nx: word count `indepv'
local nz: word count `instv'

/* numr is the number of restrictions for the F statistic */
if `"`constant'"'==`""' {
local numr `nx'-1
}
else {
local numr `nx'
}

tempname RR RM R
mat `RR'=I(`nx')
if `"`constant'"'==`""' {
mat `RM'=`RR'[1..`nx'-1,1...]
mat `R'=`RR'[1..`nx'-1,1...]
}
else {
mat `RM'=`RR'[1..`nx',1...]
mat `R'=`RR'[1..`nx',1...]
}

/* implement uhat in formular (3) in paper for kernel case */
mata: kbmscore("`depva'","`indepv'","`instv'","`touse'","`RM'",`numr') 


matrix theta_iv=thetaiv
matrix scoremm=scorem
tempname score
mat `score'=scoremm


tempname mcoef
mat `mcoef'= theta_iv


tempname thetahat
mat `thetahat'=`R'*`mcoef''


scalar signific=`level'/100
scalar nres=colsof(`mcoef')
if "`constant'" == "" {
local cv=invchi2(nres-1,signific)/*cv from chi-squared distribution with df=p, p is the number of restriction*/
local nre=colsof(`mcoef')-1
}
else {
local cv=invchi2(nres,signific)
local nre=colsof(`mcoef')
}



local sig `level'/100
/* select the alternative such that the local power is 0.75 */
local pow 0.75

/* the noncentrality parameter that gives arise power-level*/
mata:delta2(`sig',`nre',`pow')

local dlt2=delta2
local significance `level'/100
local tau 1.15
local bmin 0.00001
local bmax 1


/* implement the formular (15) in paper for kernel case, calculate the  F statistic */
mata: lrvopt("`op'",`significance',`tau',`bmin',`bmax',`cv',`dlt2',"`score'","`thetahat'") 


return matrix thetaiv=theta_iv  
return scalar F=Ft
return scalar bopt=bopt
return scalar lag=int(lg)

/* DFapprox is the second degrees of freedom */
return scalar secdf=DFapprox
return scalar adjust=adjust
}
end


mata:
void delta2(real scalar sig,
	 real scalar nre,
	 real scalar pow){
	 mm_root(x=.,&fzxp(),0,25,0.01,1000,sig,nre,pow)
	 st_numscalar("delta2",x)
	 }
     function fzxp(
	 real scalar x,      
	 real scalar sig,
	 real scalar nre,
	 real scalar pow)
{
	return(nchi2(nre,x,invchi2(nre,sig))-1+pow)
}
end


mata:
void lrvopt(string scalar op,
real scalar significance,
real scalar tau,
real scalar bmin,
real scalar bmax,
real scalar cv,
real scalar dlt2,
string scalar score,
string scalar thetahat)
{
real scalar si,tao,bi,ba,cva,delta2,gq,q,c1,c2,T,d,MBbar,a1,a2,a,bopt,Ftemp,adjust,F,DFapprox,lg
real matrix v,dep,indep,A,B,D,C,res,VA,IA,omega0,temp,omegaq1,H,Aj,MB,theta,omegahat,IROR

si=significance
tao=tau
bi=bmin
ba=bmax
cva=cv
delta2=dlt2
v=st_matrix(score)
theta=st_matrix(thetahat)

//parameter for different kernels
if (op == "B"){
gq=1
q=1
c1=1
c2=2/3
}
else if (op == "P"){
gq=6
q=2
c1=0.7500
c2=0.5393
}
else if (op == "Q"){
gq=1.4212
q=2
c1=1.2500
c2=1.0000
}

T=rows(v)
d=cols(v)


//fitting VAR(1) models
dep=v[(2::T),]
dep=dep'
indep=v[(1::T-1),]
indep=indep'
Iindep=svsolve(indep,I(rows(indep)))
A=dep*Iindep
res=dep-A*indep
VA=res*res'/(T-2*d)


IA=svsolve((I(d)-A),I(d))
// plug-in estimate of LRV
omega0=IA*VA*IA'

temp=luinv(I(rows(A))-A)
if (q==2){
omegaq1=temp*temp*temp*(A*VA+A*A*VA*A'+A*A*VA-6*A*VA*A'+VA*A'*A'+A*VA*A'*A'+VA*A')*temp'*temp'*temp'
}
else if (q==1){
H=J(rows(A),cols(A),0)
Aj=I(rows(A))
for(i=1;i<=100;i++){
H=H+Aj*VA*Aj'
Aj=Aj*A
}
omegaq1=temp*temp*A*H+H'*A'*temp'*temp'
}
// plug-in estimate of B
MB=-gq*omegaq1
// testing-optimal b
MBbar=trace(svsolve(omega0,MB))/d

if (MBbar>0){
a1=2*q*nchi2den(d,delta2,cva)*abs(MBbar)
a2=delta2*nchi2den(d+2,delta2,cva)*c2
a=a1/a2
btemp=(a^(1/(q+1)))*(T^(-q/(q+1))) //implement testing-optimal b in formular (11) in paper
}
else if (MBbar<=0){
a1=chi2den(d,cva)*cva*abs(MBbar)
a2=(tao-1)*(1-si)
a=a1/a2
btemp=(a^(1/q))/T                 //implement testing-optimal b in formular (11) in paper
}

if (btemp<=0.5){
  bopt=btemp
  }
if (btemp>=0.5){
  bopt=0.5
  }




//compute the LRV
R0=v'*v/T
R=J(d,d,0)

         
for(i=1;i<=T-1;i++){
k=i/(bopt*T)
if (op == "B" & abs(k)<=1){
R=R+(1-abs(k))*v[(i+1::T),]'*v[(1::T-i),]/T 
}
if (op == "P" & abs(k)<=0.5 & abs(k)>=0){
R=R+(1-6*k^2+6*(abs(k))^3)*v[(i+1::T),]'*v[(1::T-i),]/T 
}
else if (op == "P" & abs(k)<=1 & abs(k)>=0.5){
R=R+2*(1-abs(k))^3*v[(i+1::T),]'*v[(1::T-i),]/T 
}
if (op=="Q" & k!=0){
k=k*6*pi()/5
R=R+(3/(k^2))*(sin(k)/k-cos(k))*v[(i+1::T),]'*v[(1::T-i),]/T 
}
else if (op=="Q" & k==0){
R=R+1*v[(i+1::T),]'*v[(1::T-i),]/T 
}
}
omegahat=R0+R+R'
IROR=svsolve(omegahat,I(d))

//computer the F statistic
Ftemp=(T^(1/2))*theta'*IROR*(T^(1/2))*theta; 
adjust=0.5*(1+bopt*(c1+c2*(d-1)))+0.5*exp(bopt*(c1+c2*(d-1))) //implement the correction factor in formular (9) in paper
Ftemp=Ftemp/adjust
F=Ftemp/d
lg=bopt*T
if (op=="B"){
vdf=(floor(1/(bopt*c2)),d)
DFapprox=max(vdf)
}
else if (op=="P"|op=="Q" ){
vdf=(floor(1/(bopt*c2)),d)
DFapprox=max(vdf)-d+1
}


st_numscalar("bopt",bopt)
st_numscalar("lg",lg)
st_numscalar("Ft",F)
st_numscalar("adjust",adjust)
st_numscalar("DFapprox",DFapprox)

}
end


mata:
void kbmscore(
string scalar depva,
string scalar indepv,
string scalar instv,
string scalar touse,
string scalar RM,
real scalar numr)
{
real scalar T,m,nr
real matrix iv,iiv,Pz,indv,depv,theta_iv,res,ehat,muhat,uhat,GT,IW,IGWG


nr=numr
R=st_matrix(RM)

iv=.
st_view(iv,.,instv,touse)
indv=.
st_view(indv,.,indepv,touse)
depv=.
st_view(depv,.,depva,touse)

m=cols(iv)
T=rows(iv)
numx=cols(indv)


iiv=svsolve(iv'*iv,I(m))
Pz=iv*iiv*iv'
theta_iv=svsolve((indv'*Pz*indv),(indv'*Pz*depv))
ttheta_iv=theta_iv'
res=depv-indv*theta_iv

ehat=J(T,m,0)
for(i=1;i<=m;i++){
ehat[,i]=iv[,i]:*res
}

GT=iv'*indv/T
IW=svsolve(iv'*iv/T,I(m))
IGWG=svsolve(GT'*IW*GT,I(numx))
// uhat in formular (3) in paper
uhat=R*IGWG*GT'*IW*ehat' 
// the demean uhat
uhat=uhat'
muhat=mean(uhat,1)
muhat=muhat'
for(i=1;i<=nr;i++){
for(j=1;j<=T;j++){
uhat[j,i]=uhat[j,i]-muhat[i]
}
}

st_matrix("scorem",uhat)
st_matrix("thetaiv",ttheta_iv)

}
end


