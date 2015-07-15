/*
#################################################################################
### Jose Espinosa-Carrasco. CB/CSN-CRG. June 2015                             ### 
#################################################################################
### Code : 10.07                                                              ###
### Permutate labels of mice groups to perform different PCAs and compare     ###
### t statistic of comparisons of pairs of groups for example ts vs tseeegcg  ### 
### ./nextflow PCA_perm_test.nxf --tbl_var_MWM '*.sav'                        ###
#################################################################################
*/


params.MWM_tbl = "20150515_PCA_old_frotiersPaper/data/Ts65Dn_OLD_ACQ1_ACQ5_SUBCONJ.sav"

MWM_tbl_path = "$HOME/${params.MWM_tbl}"

println "path: $MWM_tbl_path"

//MWM_file = Channel.fromPath(MWM_tbl_path)
MWM_file = file(MWM_tbl_path)

start_perm = 1111
step = 5
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
    file MWM_file
    
    output:
    set file ('tbl_t_stat.csv') into tbl_t_stat
    set file ('tbl_t_stat_day1.csv') into tbl_t_stat_day1
    
    script:
    println "Perm is $perm"
    
    """
    export R_LIBS="/software/R/packages"
    
    Rscript \$HOME/git/mwm/lib/R/PCA_perm_test_nf.R --path2files=\$(readlink ${MWM_file}) --seed=${perm}    
    """ 
}


tbl_t_stat
    .collectFile(name: 't_stat.csv', newLine: false)
    .subscribe {
        //println "Entries are saved to file: $it"
        //println "File content is: ${it.text}"
        it.copyTo( dump_dir.resolve ( "PCA_t_statistic_${start_perm}.csv" ) )
    }
    
tbl_t_stat_day1
	.collectFile(name: 't_stat_1.csv', newLine: false)
    .subscribe {
        //println "Entries are saved to file: $it"
        //println "File content is: ${it.text}"
        it.copyTo( dump_dir.resolve ( "PCA_t_statistic_day1_${start_perm}.csv" ) )
    }
   
/*

# Rscript \$HOME/git/phecomp/lib/R/starting_regions_file_vs_24h_nf.R --tag="sum" --path2files=\$(readlink ${tr_bed_dir}) --path2plot=\$(readlink ${tr_bed_dir}) 2>&1

seed.subscribe { 
    onNext: { println it }, onComplete: { println 'Done.' } 
    }

/*
Channel
    .from( 1..3 )
    .subscribe onNext: { println it }, onComplete: { println 'Done.' }
*/
