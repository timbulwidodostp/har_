capture program drop har
program define har, eclass
	version 14,missing
	
   if !replay() {
   local vv : display "version " string(_caller()) ", missing:"
   local cmdline : copy local 0
   syntax anything [if][in],kernel(string)[/*
			*/ noCONstant  /* 
			*/ Level(cilevel) *] 
	
    _get_diopts diopts, `options'
    marksample touse	/* tvar to be subsequently marked out*/

	
    quietly{
		
	if  "`kernel'" == bsubstr("BARTLETT", 1, max(1,1)) |	///
        "`kernel'" == bsubstr("bartlett", 1, max(1,1)) |   ///
	    "`kernel'" == "bartlett" |   ///
	    "`kernel'" == "BARTLETT" |   {
        local op B
		local ktype Bartlett
	}
	
	
	if  "`kernel'" == bsubstr("PARZEN", 1, max(1,1)) |	///
	    "`kernel'" == bsubstr("parzen", 1, max(1,1)) |   ///
	    "`kernel'" == "PARZEN" |   ///
	    "`kernel'" == "parzen" |   {
		local op P
		local ktype Parzen
	}
	
	
	if  "`kernel'" == bsubstr("QUADRATIC", 1, max(1,1)) |	///
	    "`kernel'" == bsubstr("quadratic", 1, max(1,1)) |   ///
	    "`kernel'" == "QUADRATIC" |   ///
	    "`kernel'" == "quadratic" |   {
		local op Q
		local ktype Quadratic Spectral
	}
	
	if  "`kernel'" == bsubstr("ORTHOSERIES", 1, max(1,1)) |	///
	    "`kernel'" == bsubstr("orthoseries", 1, max(1,1)) |   ///
	    "`kernel'" == "ORTHOSERIES" |   ///
	    "`kernel'" == "orthoseries" |   {
		local op O
		local ktype Orthonormal Series
	}
	

	
		
			cap xt_tis
			if _rc {
				di as err "time variable not set, " /*
				*/ "use -tsset varname ...-"
				exit 111
			}
			local tvar `"`s(timevar)'"'
			markout `touse' `tvar'
			
			Checkt `tvar' `touse'
			if r(tflag)==2 {
				noi di in red `"`tvar' is regularly spaced, but does not have intervals of 1"'
				exit 198
			}
			
			if r(tflag)==1 & `"`force'"' == "" {
				noi di in red `"`tvar' is not regularly spaced"'
				exit 198
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
			
	
	
	_iv_parse `0'


	local lhs `s(lhs)'
	local endogc `s(endog)'
	local exogc `s(exog)'
	local instc `s(inst)'

	
	tempvar normwt 
	qui gen double `normwt' = 1 if `touse'

	
	/* Remove collinear variables */
	CheckCollin `lhs' if `touse' [iw=`normwt'],endog(`endogc') exog(`exogc') inst(`instc') `perfect'
	
	local endog `s(endog)'
	local exog `s(exog)'
	local inst `s(inst)'
	
		NoOmit `endog', touse(`touse')
		tempname Dendog
		matrix `Dendog'=r(noomit)
		
		NoOmit `exog', touse(`touse')
		tempname Dexog
		matrix `Dexog'=r(noomit)
		
		NoOmit `inst', touse(`touse')
		tempname Dinst
		matrix `Dinst'=r(noomit)
		
		scalar nendo=0
		local newendog `""'
		foreach var of local endog {
		scalar nendo=nendo+1
		if (`Dendog'[1,nendo]!=0){
		local newendog  `"`newendog' `var'"'
		}
		else{
		local newendog  `"`newendog'"'
		}
		}
		
		
		scalar nexog=0
		local newexog `""'
		foreach var of local exog {
		scalar nexog=nexog+1
		if (`Dexog'[1,nexog]!=0){
		local newexog  `"`newexog' `var'"'
		}
		else{
		local newexog  `"`newexog'"'
		}
		}
		
		
		scalar ninst=0
		local newinst `""'
		foreach var of local inst {
		scalar ninst=ninst+1
		if (`Dinst'[1,ninst]!=0){
		local newinst  `"`newinst' `var'"'
		}
		else{
		local newinst  `"`newinst'"'
		}
		}
		
	
	local endog `newendog'
	local exog  `newexog'
	local inst  `newinst'

    local full `exogc' `endogc' `instc'
	local sub  `exog' `endog' `inst'
	local remove: list full -sub
	local nremove: word count `remove'	
	
    local depva `lhs'
	local indepv `"`endog' `exog' `CONS'"'
	local instv `"`inst' `exog' `CONS'"'
	local nx:word count `indepv'
	local nz:word count `instv'
	local nendog: word count `endog'
	
	
	

	
			
	    if("`op'"=="O"){
	
	    if(`nendog'!=0){
        ivregress gmm `depva' `exog' (`endog'=`inst')  if `touse',`carg' 
	    }
        else {
        ivregress gmm `depva' `exog'   if `touse',`carg' 
	    }
			if e(N)==0 | e(N)>=. { 
				di in red `"no observations"'
				exit 2000
            }
			
	 		local nobs=e(N)
			local mdf=e(df_m)
			local tdf=e(df_r)
			
			tempname beta 
			mat `beta' = e(b)	
		
			
		 /* compute F statistic for OS case */ 
		 if(`nendog'!=0){
         noi OS_F_all `depva' `exog' (`endog'=`inst')  if `touse',level(`level') `carg' 
	     }
         else {
         noi OS_F_all `depva' `exog'   if `touse',level(`level') `carg' 
	     }
		  

			tempname betahat
			mat `betahat'=r(thetaiv) /* the IV coefficient vector */
			scalar msecdf=r(secdf)   /* second degrees of freedom for F statistic */
			scalar mF=r(F)           /* the F statistic */
			scalar kopt=r(kopt)      /* the data-driven optimal K */
			
						
	
		/* compute t statistic for OS case */ 	
		if(`nendog'!=0){
        noi OS_t `depva' `exog' (`endog'=`inst')  if `touse',level(`level') `carg' 
	    }
        else {
        noi OS_t `depva' `exog'   if `touse',level(`level') `carg' 
	    }
			tempname  udf ut uadjstderror ubetahat
			mat `ubetahat'=r(uthetaiv)        /* the IV coefficient vector */
			mat `uadjstderror'=r(adjstderror) /* the adjusted std error    */
			mat `udf'=r(udf)                  /* the degrees of freedom for t statistic */
			mat `ut'=r(ut)                    /* the t statistic */                    
 	  
	  }
			 		 
	   
	   if ("`op'"=="B" |"`op'"=="P"|"`op'"=="Q"){

	   if(`nendog'!=0){
       ivregress gmm `depva' `exog' (`endog'=`inst')  if `touse',`carg' 
	   }
       else {
       ivregress gmm `depva' `exog'   if `touse',`carg' 
	   }
			
			if e(N)==0 | e(N)>=. { 
				di in red `"no observations"'
				exit 2000
			}
			local nobs=e(N)
			local mdf=e(df_m)
			local tdf=e(df_r)
						


			tempname beta 
			mat `beta' = e(b)
	  
	     /* compute the F statistic for kernel case */
		 if(`nendog'!=0){
         noi kernel_F_all `depva' `exog' (`endog'=`inst')  if `touse',kernel(`op') level(`level') `carg'
	     }
         else {
         noi kernel_F_all `depva' `exog'   if `touse',kernel(`op') level(`level') `carg'
	     }

			tempname kbetahat
			mat `kbetahat'=r(thetaiv)   /* the IV coefficient vector */
			scalar kmbopt=r(bopt)       /* the data-driven b */
			scalar kmF=r(F)             /* the F statistic */
			scalar kmsecdf=r(secdf)     /* second degrees of freedom for F statistic */
			scalar kmlag=r(lag)         /* the data-driven truncation lag */
		
		
		
		
		/* compute the t statistic for kernel case */
		 if(`nendog'!=0){
         noi kernel_t `depva' `exog' (`endog'=`inst')  if `touse',kernel(`op') level(`level') `carg'
	     }
         else {
         noi kernel_t `depva' `exog'   if `touse',kernel(`op') level(`level') `carg'
	     }
			
			tempname  kut kudf  kuadjstderror kubetahat
			mat `kubetahat'=r(uthetaiv)          /* the IV coefficient vector */
            mat `kut'=r(ut)                      /* the t statistic */
            mat `kudf'=r(udf)                    /* the degrees of freedom for t statistic */
			mat `kuadjstderror'=r(adjstderror)   /* the adjusted std error */
					
			}
			
			
		if(`nendog'!=0){
        local 0 `depva' `exog' (`endog'=`inst') 
	    }
        else {
        local 0 `depva' `exog' 
	    }
			
          /* the first degrees of freedom of F statistic for all kernels and orthonormal series */
		  eret scalar fdf=`mdf'
		  
	      if ("`op'"=="O" ){
		    eret scalar sF=mF
			eret scalar ssdf=msecdf
			eret matrix sbetahat=`betahat'
			eret matrix st=`ut'
			eret matrix sdf=`udf'
			eret matrix sstderr=`uadjstderror'
			eret scalar kopt=kopt
			eret local depvar `"`depva'"'
			eret scalar N=`nobs'
			eret local vcetype `"HAR"'
			eret local kerneltype `"`ktype'"'
			
			

			version 10: ereturn local cmdline `"har `cmdline'"'
			local varline `0'
			eret local varline  `varline'
            eret local carg    `carg'
			eret local cmd     `"har"'
			eret local title /*
				*/ "Regression with HAR standard errors"
			global S_E_var `"`e(varline)'"'
			_post_vce_rank
          
		   }
		   
		
		
		else{ 
		
			eret scalar kF=kmF
			eret scalar ksdf=kmsecdf
			eret matrix kbetahat=`kbetahat'
			eret matrix kt=`kut'
			eret matrix kdf=`kudf'
			eret matrix kstderr=`kuadjstderror'
	   		eret local  depvar `"`depva'"'
			eret scalar N =  `nobs'
	        eret scalar lag =kmlag    
		    eret local  vcetype `"HAR"'
			eret local  kerneltype `"`ktype'"'
			


			version 10: ereturn local cmdline `"har `cmdline'"'
			local varline `0'
			eret local varline  `varline'
			eret local carg    `carg'
			eret local cmd     `"har"'
			eret local title /*
				*/ "Regression with HAR standard errors"
			global S_E_var `"`e(varline)'"'
			_post_vce_rank
         
		 }     
    	 }
	     }
		else {
	    if `"`e(cmd)'"' != `"har"' {
			error 301 
		}
		syntax,kernel(string) [, Level(cilevel) *]
		_get_diopts diopts, `options'
	   }
	
	 if ("`op'"=="O" ){
	 	 if (`nremove'!=0){
		 di "{txt}note: `remove' dropped because of collinearity"
		 }
		 #delimit;
		 di in gr `"Regression with HAR standard errors"'
		_col(53)
		`"Number of obs  ="' in yel %10.0f e(N) _n
		 in gr `"Kernel: "' in ye e(kerneltype) 
	   	_col(53) 
		in gr 
		`"F("' in gr %3.0f e(fdf) in gr `","' in gr %6.0f e(ssdf)
		in gr `")"' _col(68) `"="' in ye %10.2f e(sF) _n
		in gr `"Data-driven optimal K: "' in ye e(kopt)
		_col(53) in gr `"Prob > F       =    "' 
		in ye %6.4f fprob(e(fdf),e(ssdf),e(sF)) _n ;
		#delimit cr 
		}
		 else{
		 if (`nremove'!=0){
		 di "{txt}note: `remove' dropped because of collinearity" 
		 }
		 #delimit;
		 di in gr `"Regression with HAR standard errors"'
		_col(53)
		`"Number of obs  ="' in yel %10.0f e(N) _n
		 in gr `"Kernel: "' in ye e(kerneltype) 
	     _col(53) 
		in gr 
		`"F("' in gr %3.0f e(fdf) in gr `","' in gr %6.0f e(ksdf)
		in gr `")"' _col(68) `"="' in ye %10.2f e(kF) _n
		in gr `"Data-driven optimal lag: "' in ye e(lag)
		_col(53) in gr `"Prob > F       =    "' 
		in ye %6.4f fprob(e(fdf),e(ksdf),e(kF)) _n ;
	    #delimit cr
		}	 


   
   
   har_tab,kernel("`kernel'") level(`level') `diopts'


	
end

capture program drop Checkt
program define Checkt, rclass sort
	args tvar touse 

	replace `touse'=. if `touse'==0
	summ `touse', meanonly
	if r(N) == 0 {
		exit 2000
	}
	ret scalar tflag = 0 
	sort `touse' `tvar'
	tempname tsdelta
	if "`: char _dta[_TSdelta]'" != "" {
		scalar `tsdelta' = `: char _dta[_TSdelta]'
	}
	else {
		scalar `tsdelta' = 1
	}
	tempvar tt
	gen double `tt' = (`tvar'-`tvar'[_n-1])/`tsdelta' if `touse'<.	
	summ `tt', meanonly
	if r(min) != r(max) {
		ret scalar tflag = 1
	}
	if r(mean) > 1 & r(min) == r(max) {
		ret scalar tflag = 2
	}
	replace `touse'=0 if `touse'>=.
end

capture program drop CheckCollin
program CheckCollin, sclass
	if _caller() >= 11 {
		local vv : di "version " string(_caller()) ":"
	}
	syntax varlist(ts min=1 max=1) [if] [in] [iw/]	 	///
		[, endog(varlist fv ts) exog(varlist fv ts) 	///
		   inst(varlist fv ts) PERfect NOCONSTANT ]
	marksample touse
	if `"`exp'"' != "" {
		local wgt `"[`weight'=`exp']"'
	}
        local fvops = "`s(fvops)'" == "true" | _caller() >= 11
        local tsops = "`s(tsops)'" == "true" 

        if `fvops' {
                if _caller() < 11 {
                        local vv "version 11:"
                }
		local expand "expand"
		fvexpand `exog' if `touse'
		local exog  "`r(varlist)'"
		fvexpand `inst' if `touse'
		local inst  "`r(varlist)'"
		fvexpand `endog'
		local endog "`r(varlist)'"
	}
	/* Catch specification errors */	
	/* If x in both exog and endog, error out */
	local both : list exog & endog
	foreach x of local both {
		di as err 	///
"`x' included in both exogenous and endogenous variable lists"
		exit 498
	}
	
	if "`perfect'" == "" {
		/* If x in both endog and inst, error out */
		local both : list endog & inst
		foreach x of local both {
			di as err 	///
"`x' included in both endogenous and excluded exogenous variable lists"
			exit 498
		}
	}
	
	/* If x on both LHS and (RHS or inst), error out */
	local both : list varlist & endog
	if "`both'" != "" {
		di as err 	///
		 "`both' specified as both regressand and endogenous regressor"
		exit 498
	}
	local both : list varlist & exog
	if "`both'" != "" {
		di as err 	///
		   "`both' specified as both regressand and exogenous regressor"
		exit 498
	}

	local both : list varlist & inst
	if "`both'" != "" {
		di as err 	///
"`both' specified as both regressand and excluded exogenous variable"
		exit 498
	}

	/* Now check for collinearities */
	`vv' ///
	_rmdcoll `varlist' `endog' `exog' `wgt' if `touse', `noconstant'
	local totvarlist  `r(varlist)'
	if "`r(k_omitted)'" == "" {
		local both `r(varlist)'
		local endog : list endog & both
		local exog  : list exog & both
	}
	else {
		local list `r(varlist)'
		local omitted `r(k_omitted)'
		if `omitted' {
			foreach var of local list {
				_ms_parse_parts `var'
				local inendog : list var in endog
				local inexog : list var in exog
				if (`inendog') {
					local endog_keep `endog_keep' `var'
				}
				else {
					local exog_keep `exog_keep' `var'
				}
			}
		}
		else {
			local exog_keep `exog'
			local endog_keep `endog'
		}
		local endog `endog_keep'
		local exog `exog_keep'
	}
	`vv' ///
	_rmcoll `inst', `expand' `noconstant'
	local inst `r(varlist)'
	if "`inst'" != "" & "`endog'`exog'" != "" {
		if "`noconstant'"  == "" {
			tempvar tmpcons
			qui gen double `tmpcons' = 1 if `touse'
		}
		else {
			local tmpcons
		}
		if "`perfect'" == "" {
			`vv' ///
			_rmcoll2list, alist(`endog' `exog' `tmpcons') ///
				blist(`inst') ///
				normwt(`exp') touse(`touse')
			local inst `r(blist)'
		}
		else if "`exog'" != "" {   // allowing perfect instruments
			`vv' ///
			_rmcoll2list, alist(`exog' `tmpcons') blist(`inst') ///
				normwt(`exp') touse(`touse')
			local inst `r(blist)'
		}
	}
	sreturn local endog `endog'
	sreturn local exog `exog'
	sreturn local inst `inst'
	sreturn local fvops `fvops'
	sreturn local tsops `tsops'
	
end

capture program drop NoOmit
program NoOmit, rclass
	syntax [varlist(fv ts default=none)] [,touse(string) exporder(string)]
	if "`varlist'" == "" {
		local hasfv = 0
		local omitted 0 
		local noomitted 0
		return scalar omitted = `omitted'
		return scalar noomitted = `noomitted'
		return scalar hasfv = `hasfv'
		exit
	}
	local hasfv = "`s(fvops)'" == "true"	
	if ("`exporder'" == "") {
		fvexpand `varlist' if `touse'
		local full_list `r(varlist)'
	}
	else {
		local full_list `exporder'
	}
	local cols : word count `full_list'
	tempname noomit noomitcols
	mat `noomit' = J(1,`cols',1)
	local omitted 0
	local noomitted 0
	local i 1
	foreach var of local full_list {
		_ms_parse_parts `var'
		if `r(omit)' {
			local ++omitted	
			mat `noomit'[1,`i'] == 0
			if !`hasfv' {
				local hasfv 1
			}
		}
		else {
			local ++noomitted
			capture assert `noomitcols'[1,1]>0
			if !_rc {
				mat `noomitcols' = `noomitcols'[1,1...],`i'
			}
			else {
				mat `noomitcols' = J(1,1,`i')
			}
		}
		local ++i
	}
	return scalar omitted = `omitted'
	return scalar noomitted = `noomitted'
	if `noomitted' > 0 {
		return matrix noomitcols = `noomitcols'
	}
	return scalar hasfv = `hasfv'
	return matrix noomit = `noomit'
end





