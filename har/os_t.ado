program define OS_t, rclass
/*define the syntax*/
syntax anything [if][in],[,noConstant Level(cilevel)]
marksample touse
qui count if `touse'
scalar numobs=r(N)

qui{

	_iv_parse `0'

	local lhs `s(lhs)'
	local endog `s(endog)'
	local exog `s(exog)'
	local inst `s(inst)'
    local 0 `s(zero)'
	
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
	local instv `inst' `exog' `CONS'

	
	
local nx: word count `indepv'
local nz: word count `instv'

tempname RR RM R
mat `RR'=I(`nx')




mat udf=J(1,`nx',0)
mat uF=J(1,`nx',0)
mat uadjstderror=J(1,`nx',0)
mat ut=J(1,`nx',0)


forvalues i = 1/`nx'{
mat `R'=`RR'[`i'..`i',1...]
mat `RM'=`RR'[`i'..`i',1...]

local numr 1


/* implement uhat in formular (3) in paper */
mata: umscore("`depva'","`indepv'","`instv'","`touse'","`RM'",`numr') 

matrix theta_iv=thetaiv
matrix scoremm=scorem
tempname score
mat `score'=scoremm

tempname mcoef
mat `mcoef'= theta_iv


tempname thetahat
mat `thetahat'=`R'*`mcoef''

scalar signific=`level'/100
local cv=invchi2(1,signific)



local nre=1
local sig `level'/100
/* select the alternative such that the local power is 0.75 */
local pow 0.75
/* the noncentrality parameter that gives arise power-level*/
mata: sudelta2(`sig',`nre',`pow')


local dlt2=delta2
local tau 1.15
local significance `level'/100

/* calculate the F(t_square) statistic */
mata: ulrvseries(`nre',`significance',`tau',`cv',`dlt2',"`score'","`thetahat'") 

/* for univariate t test, the degrees of freedom is kpot  */
mat udf[1,`i']=ukopt
mat uF[1,`i']=F
mat uadjstderror[1,`i']=sqrt(omegahat)
mat ut[1,`i']= sign(theta_iv[1,`i'])*sqrt(F)
}

return matrix uthetaiv=theta_iv
return matrix udf=udf
return matrix adjstderror=uadjstderror
return matrix ut=ut
}
end


mata:
void ulrvseries(
real scalar nre,
real scalar significance,
real scalar tau,
real scalar cv,
real scalar dlt2,
string scalar score,
string scalar thetahat)
{
real scalar T,d,MBbar,nr,kstar,cns,cva,tao,delta2,sig,Ftemp,F,k
real matrix v,dep,indep,A,B,D,C,res,VA,IA,omega0,temp,omegaq1,MB,dftv,omegahat,theta

nr=nre
sig=significance
tao=tau
cva=cv
delta2=dlt2
v=st_matrix(score)
theta=st_matrix(thetahat)

T=rows(v)
d=cols(v)


// fast fourier transform
dftv=J(T,d,1i)
for(i=1;i<=d;i++){
x=v[(1::T),i]
y=C(x)
yr=Re(y)
yi=Im(y)
n=rows(x)
dftr=J(n,1,.)
dfti=J(n,1,.)
range=range(1,n,1)
for (j=0; j<=n-1; j++) {	
ae=-2*pi()/n
dftr[j+1]=colsum(yr:*cos((range:-1):*ae:*j)-yi:*sin((range:-1):*ae:*j))
dfti[j+1]=colsum(yr:*sin((range:-1):*ae:*j)+yi:*cos((range:-1):*ae:*j))
}
dftv[(1::T),i]=conj(C(dftr,dfti))
}


//VAR(1) plug in
dep=v[(2::T),]
dep=dep'
indep=v[(1::T-1),]
indep=indep'
Iindep=svsolve(indep,I(rows(indep)))
A=dep*Iindep
res=dep-A*indep
VA=res*res'/(T-nr)


IA=svsolve((I(d)-A),I(d))

// plug-in estimate of the LRV
omega0=IA*VA*IA'
temp=luinv(I(rows(A))-A)
omegaq1=temp*temp*temp*(A*VA+A*A*VA*A'+A*A*VA-6*A*VA*A'+VA*A'*A'+A*VA*A'*A'+VA*A')*temp'*temp'*temp'

// plug-in estimate of B
MB=-(pi()^2)/6*omegaq1


//testing-optimal K
MBbar=trace(MB*svsolve(omega0,I(d)))/d
if (MBbar>0){
a1=4*nchi2den(d,delta2,cva)*abs(MBbar)
a2=delta2*nchi2den(d+2,delta2,cva)*1
a=a2/a1
Ktemp=(a^(1/3))*(T^(2/3))
}
else if (MBbar<=0){
a1=chi2den(d,cva)*cva*abs(MBbar)
a2=(tao-1)*(1-sig)
a=a2/a1
Ktemp=(a^(1/2))*T
}

 if (Ktemp<=nr+4) {
  Kwtemp=nr+4
 }
 else if (Ktemp>nr+4 & Ktemp<=T){
  Kwtemp=Ktemp
 }
 else if(Ktemp>T){
  Kwtemp=T
 }
  
   maxkstar=(floor(Kwtemp/2),1)
   kstar=max(maxkstar)
   k=2*kstar
   
  

omegahat=Re(dftv[2::kstar+1,]'*dftv[2::kstar+1,]/(kstar*T))



//compute the F statistic
Ftemp=sqrt(T)*theta'*svsolve(omegahat,I(d))*sqrt(T)*theta

//adjusted F statisic
F=((k-nr+1)/(nr*k))*Ftemp
omegahat=omegahat/T*((nr*k)/(k-nr+1))
st_numscalar("omegahat",omegahat)
st_numscalar("ukopt",k)
st_numscalar("F",F)
}
end

mata:
void sudelta2(real scalar sig,
	 real scalar nre,
	 real scalar pow){
	 mm_root(x=.,&sufzxp(),0,25,0.01,1000,sig,nre,pow)
	 st_numscalar("delta2",x)
	 }
     function sufzxp(
	 real scalar x,      
	 real scalar sig,
	 real scalar nre,
	 real scalar pow)
{
	return(nchi2(nre,x,invchi2(nre,sig))-1+pow)
}
end


mata:
void umscore(
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



