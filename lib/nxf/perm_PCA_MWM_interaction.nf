/*
#################################################################################
### Jose Espinosa-Carrasco. CB/CSN-CRG. March 2016                            ### 
#################################################################################
### Code : 14.03                                                              ###
### Permutate labels of mice groups to perform different PCAs and compare     ###
### t statistic of comparisons of pairs of interactions of group and day      ### 
### for example ts_day4 vs ts_day5 or ts_eeegcg_day4 vs ts_day5               ###
### Young animals acquisition                                                 ### 
### ./nextflow perm_PCA_MWM_interaction.nxf                                   ###
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
    set file ('tbl_t_stat_interaction.csv') into tbl_t_stat
    //set file ('tbl_t_stat_rem1.csv') into tbl_t_stat_day1
    
    script:
    println "Perm is $perm"
    
    """
    export R_LIBS="/software/R/packages"
    
    Rscript \$HOME/git/mwm/lib/R/PCA_perm_test_nf_interactionGroupDay.R --path2files=\$(readlink ${MWM_file}) --seed=${perm}    
    """ 
}


tbl_t_stat
    .collectFile(name: 't_stat_interaction.csv', newLine: false)
    .subscribe {
        //println "Entries are saved to file: $it"
        //println "File content is: ${it.text}"
        it.copyTo( dump_dir.resolve ( "PCA_t_statistic_interaction_${start_perm}.csv" ) )
    }