* har_tab.ado
capture program drop har_tab
program har_tab, eclass
         syntax,kernel(string)[, Level(cilevel)]
		 tempname mytab z b t p ll ul adjsterr df
        .`mytab' = ._tab.new, col(8) lmargin(0)
        .`mytab'.width     10    |11   11     8    6     8     12    12
        .`mytab'.titlefmt  .     .     .     %6s   .     .     %24s  .
        .`mytab'.pad       .     2     1      0    0     2     3     3
        .`mytab'.numfmt    .    %9.0g %9.0g  %7.2f %5.0f %5.3f %9.0g %9.0g

        local namelist : colname e(b)
        local eqlist : coleq e(b)
        local k : word count `namelist'
        .`mytab'.sep, top
        if `:word count `e(depvar)'' == 1 {
                local depvar "`e(depvar)'"
        }
        
 		
		.`mytab'.titles ""                      /// 1
                        ""                      /// 2
	                   "HAR   "                    /// 3
                        ""                      /// 4
						""                      /// 5
                        ""                      /// 6
                        "" ""                   //  7 8
	
	
		.`mytab'.titles "`depvar'"                      /// 1
                        "Coef."                         /// 2
	                   "Std.Err."                       /// 3
                        "t"                             /// 4
						"df"                            /// 5
                        "P>|t|"                         /// 6
                        "[`level'% Conf. Interval]" ""  //  7 8
    	
		
		forvalues i = 1/`k' {
                local name : word `i' of `namelist'
                local eq   : word `i' of `eqlist'
                if "`eq'" != "_" {
                        if "`eq'" != "`eq0'" {
                                .`mytab'.sep
                                  local eq0 `"`eq'"'
                                .`mytab'.strcolor result  .  .  .  .  .  .  .
                                .`mytab'.strfmt    %-12s  .  .  .  .  .  .  .
                                .`mytab'.row      "`eq'" "" "" "" "" "" "" ""
                                .`mytab'.strcolor   text  .  .  .  .  .  .  .
                                .`mytab'.strfmt     %12s  .  .  .  .  .  .  .
                        }
                        local beq "[`eq']"
						                }
                else if `i' == 1 {
                        local eq
                        .`mytab'.sep
                }

				
	   if ("`kernel'"=="orthoseries"| "`kernel'"=="ORTHOSERIES" | "`kernel'"=="O" |"`kernel'"=="o"){
		  
		  
                scalar `t' = el(e(st),1,`i')
                local tstat = `t'
				scalar `df'=el(e(sdf),1,`i')
				scalar `b'=el(e(sbetahat),1,`i')
				scalar `adjsterr'=el(e(sstderr),1,`i')
	            scalar `p' = 2*ttail(`df',(abs(`tstat')))
                scalar `ll' = `b'-`adjsterr'*invttail(`df',((100-`level')/200))
                scalar `ul' = `b'+`adjsterr'*invttail(`df',((100-`level')/200))
                .`mytab'.row    "`name'"                ///
                                `b'                    ///
                                `adjsterr'          ///
                                `t'                 ///
                                `df'                ///
							    `p'                     ///
                                `ll' `ul'
        }
				
		else {
			
			
                scalar `t' =el(e(kt),1,`i')
                local tstat = `t'
				scalar `df'=el(e(kdf),1,`i')
				scalar `b'=el(e(kbetahat),1,`i')
				scalar `adjsterr'=el(e(kstderr),1,`i')
	            scalar `p' = 2*ttail(`df',(abs(`tstat')))
                scalar `ll' = `b'-`adjsterr'*invttail(`df',((100-`level')/200))
                scalar `ul' = `b'+`adjsterr'*invttail(`df',((100-`level')/200))
                .`mytab'.row    "`name'"                ///
                                `b'                   ///
                                `adjsterr'          ///
                                `t'                 ///
                                `df'                ///
							    `p'                     ///
                                `ll' `ul'
		 
      }
		 
	  }
        .`mytab'.sep, bottom
		

end

