/*
#################################################################################
### Jose Espinosa-Carrasco. CB/CSN-CRG. March 2016                            ### 
#################################################################################
### Code : 14.03                                                              ###
### Permutate labels of mice groups to perform different PCAs and compare     ###
### t statistic to compare PC2 of treated vs untreated genotypes on           ### 
### acquistion day 4 and 5                                                    ###
### Young animals acquisition                                                 ### 
### ./nextflow perm_PCA_MWM_PC2_day4_5.nf                                     ###
#################################################################################
*/

params.MWM_tbl = "20151001_ts65_young_MWM/data/ts65_young.csv"
MWM_tbl_path = "$HOME/${params.MWM_tbl}"

println "path: $MWM_tbl_path"

MWM_file = file(MWM_tbl_path)

start_perm = 1111
step = 10000
end_perm = start_perm + step
perm = Channel.from(start_perm..end_perm)

dump_dir = file("$HOME/20151001_ts65_young_MWM/data/")

process perm {
    input:
    val perm from perm
    file MWM_file
    
    output:
    set file ('tbl_t_stat_PC2_day4_vs_5.csv') into tbl_t_stat
        
    script:
    println "Perm is $perm"
    
    """
    export R_LIBS="/software/R/packages"
    
    Rscript \$HOME/git/mwm/lib/R/PCA_perm_test_day_4_5_PC2_young_nf.R --path2files=\$(readlink ${MWM_file}) --seed=${perm}    
    """ 
}


tbl_t_stat
    .collectFile(name: 't_stat_interaction.csv', newLine: false)
    .subscribe {
        //println "Entries are saved to file: $it"
        //println "File content is: ${it.text}"
        it.copyTo( dump_dir.resolve ( "PCA_t_statistic_PC2_day4_5_${start_perm}.csv" ) )
    }