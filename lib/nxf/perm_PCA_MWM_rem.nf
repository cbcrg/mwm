/*
#################################################################################
### Jose Espinosa-Carrasco. CB/CSN-CRG. July 2015                             ### 
#################################################################################
### Code : 10.07                                                              ###
### Permutate labels of mice groups to perform different PCAs and compare     ###
### t statistic of comparisons of pairs of groups for example ts vs tseeegcg  ###
### REMOVAL SESSION                                                           ### 
### ./nextflow PCA_perm_test_rem.nxf                                          ###
#################################################################################
*/


start_perm = 1111
step = 10

end_perm = start_perm + step
perm = Channel.from(start_perm..end_perm)

dump_dir = file("$HOME/git/mwm/lib/nxf/")
//seed = Channel.from(111,222,333)

/*
perm
    .subscribe {
        onNext: { println it } 
        onComplete: { println 'Done.' } 
    }
*/
      
process perm {
    input:
    val perm from perm
    
    output:
    set file ('tbl_t_stat.csv') into tbl_t_stat
    
    script:
    println "Perm is $perm"
    
    """
    export R_LIBS="/software/R/packages"
    
    Rscript \$HOME/git/mwm/lib/R/PCA_perm_test_rem_nf.R --seed=${perm}    
    """ 
}


tbl_t_stat
    .collectFile(name: 't_stat.csv', newLine: false)
    .subscribe {
        //println "Entries are saved to file: $it"
        //println "File content is: ${it.text}"
        it.copyTo( dump_dir.resolve ( "PCA_t_statistic_${start_perm}.csv" ) )
    }
