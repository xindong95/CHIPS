#MODULE: bin based methods

#PARAMETERS
# _logfile=output_path + "/logs/peaks.log"
_bin_size="100"

# # for narrowpeak calling, extra parameters passed to macs2.
# # e.g, --nomodel.  --nomodel is turned on for broad peaks by default,
# # for narrowpeak, has to modify the below parameter.
# _macs_extra_param = config.get("macs_extra_param", "")

def getTreats(wildcards):
    tmp=[]
    #print(wildcards)
    rep_s = wildcards.rep
    rep_n = int(rep_s[-1])
    #print(rep_n)

    #USE the formula: treatment = 2^(rep-1); control = treatment+1
    treat = 2**(rep_n-1) if rep_n > 1 else 0
    r = config['runs'][wildcards.run]
    #print(r)
    #print(treat)
    if treat < len(r) and r[treat]:
        tmp = [output_path + "/align/%s/%s_unique.sorted.bam" % (r[treat],r[treat])]
    #print("TREAT: %s" % tmp)
    if not tmp:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! Missing treatment file for run: %s, rep: %s" % (wildcards.run, rep_s)]
    return tmp

def getConts(wildcards):
    tmp=[]
    #print(wildcards)
    rep_s = wildcards.rep
    rep_n = int(rep_s[-1])

    #USE the formula: treatment = 2^(rep-1); control = treatment+1
    cont = 2**(rep_n-1) + 1 if rep_n > 1 else 1
    r = config['runs'][wildcards.run]
    #print(r)
    #print(cont)
    if cont < len(r) and r[cont]:
        tmp = [output_path + "/align/%s/%s_unique.sorted.bam" % (r[cont],r[cont])]
    #print("CONT: %s" % tmp)
    return tmp

def getFilteredTreats(wildcards):
    tmp=[]
    #print(wildcards)
    rep_s = wildcards.rep
    rep_n = int(rep_s[-1])
    #print(rep_n)

    #USE the formula: treatment = 2^(rep-1); control = treatment+1
    treat = 2**(rep_n-1) if rep_n > 1 else 0
    r = config['runs'][wildcards.run]
    #print(r)
    #print(treat)
    if r[treat]:
        treatSample = config["samples"][r[treat]]
        if len(treatSample) > 1 and ('cutoff' in config) and config['cutoff']:
            #PE
            if treat < len(r):
                tmp = [output_path + "/align/%s/%s_unique.sorted.sub%s.bam" % (r[treat],r[treat],config['cutoff'])]
            #print("TREAT: %s" % tmp)
            if not tmp:
                #NOTE: I can't figure out a proper kill command so I'll do this
                tmp=["ERROR! Missing treatment file for run: %s, rep: %s" % (wildcards.run, rep_s)]
        # else:
        #     #SE
        #     if treat < len(r):
        #         tmp = [output_path + "/align/%s/%s_unique.sorted.bam" % (r[treat],r[treat])]
            #print("TREAT: %s" % tmp)
            # if not tmp:
            #     #NOTE: I can't figure out a proper kill command so I'll do this
            #     tmp=["ERROR! Missing treatment file for run: %s, rep: %s" % (wildcards.run, rep_s)]
    return tmp

def getFilteredConts(wildcards):
    tmp=[]
    #print(wildcards)
    rep_s = wildcards.rep
    rep_n = int(rep_s[-1])

    #USE the formula: treatment = 2^(rep-1); control = treatment+1
    cont = 2**(rep_n-1) + 1 if rep_n > 1 else 1
    r = config['runs'][wildcards.run]
    #print(r)
    #print(cont)
    if r[cont]:
        contSample = config["samples"][r[cont]]
        if len(contSample) > 1 and ('cutoff' in config) and config['cutoff']:
            #PE
            if cont < len(r):
                tmp = [output_path + "/align/%s/%s_unique.sorted.sub%s.bam" % (r[cont],r[cont],config['cutoff'])]
        # else:
        #     #SE
        #     if cont < len(r):
        #         tmp = [output_path + "/align/%s/%s_unique.sorted.bam" % (r[cont],r[cont])]
        # #print("CONT: %s" % tmp)
    return tmp

def checkBAMPE(wildcards):
    """Fn returns '-f BAMPE' IF the run's FIRST treatment replicate (sample) is
    Paired-END.
    NOTE: this flag is required for macs2 callpeak, AUTO detect format does not
    work with PE bams!
    """
    r = config['runs'][wildcards.run]
    #GET first treatement sample
    first = config['samples'][r[0]]
    ret = "-f BAMPE" if len(first) == 2 else ""
    return ret

#NOTE: using the _refs from chips.snakefile
def bin_targets(wildcards):
    """Generates the targets for this module"""
    #print(wildcards)
    ls = []
    for run in config["runs"].keys():
        #NOTE: using the fact that _reps has this info parsed already!
        for rep in _reps[run]:
            #GENERATE Run name: concat the run and rep name
            runRep = "%s.%s" % (run, rep)
    #         if ("macs2_broadpeaks" in config) and config["macs2_broadpeaks"]:
    #             ls.append(output_path + "/peaks/%s/%s_sorted_peaks.broadPeak" % (runRep,runRep))
    #             ls.append(output_path + "/peaks/%s/%s_sorted_peaks.bed" % (runRep,runRep))
    #             # ls.append(output_path + "/peaks/%s/%s_sorted_summits.bed" % (runRep,runRep))
    #         else:
    #             ls.append(output_path + "/peaks/%s/%s_sorted_peaks.narrowPeak" % (runRep,runRep))
    #             ls.append(output_path + "/peaks/%s/%s_sorted_peaks.bed" % (runRep,runRep))
    #             ls.append(output_path + "/peaks/%s/%s_sorted_summits.bed" % (runRep,runRep))
    #         ls.append(output_path + "/peaks/%s/%s_5fold_peaks.bed" % (runRep,runRep))
    #         ls.append(output_path + "/peaks/%s/%s_peaks.png" % (runRep,runRep))
    #         ls.append(output_path + "/peaks/%s/%s_treat_pileup.bw" % (runRep,runRep))
    #         ls.append(output_path + "/peaks/%s/%s_control_lambda.bw" % (runRep,runRep))
    #         ls.append(output_path + "/peaks/%s/%s_treat_pileup.sorted.bdg.gz" % (runRep,runRep))
    #         ls.append(output_path + "/peaks/%s/%s_control_lambda.sorted.bdg.gz" % (runRep,runRep))
    #         ls.append(output_path + "/peaks/%s/%s_treatment.igv.xml" % (runRep,runRep))
    #         ls.append(output_path + "/peaks/%s/%s_peaks.bed" % (runRep,runRep))
    #         ls.append(output_path + "/peaks/%s/%s_model.R" % (runRep,runRep))
    #         if ('cutoff' in config) and config['cutoff'] and len(config['samples'][config['runs'][run][0]]) == 2:
    #             if ("macs2_broadpeaks" not in config) or (config["macs2_broadpeaks"] == False):
    #                 ls.append(output_path + "/peaks/%s/%s.sub%s_summits.bed" % (runRep,runRep,config["cutoff"]))
    #             ls.append(output_path + "/peaks/%s/%s.sub%s_peaks.bed" % (runRep,runRep,config["cutoff"]))
    #             ls.append(output_path + "/peaks/%s/%s.sub%s_treat_pileup.bw" % (runRep,runRep,config["cutoff"]))
    #             ls.append(output_path + "/peaks/%s/%s.sub%s_control_lambda.bw" % (runRep,runRep,config["cutoff"]))
    # ls.append(output_path + "/peaks/peakStats.csv")
    # ls.append(output_path + "/peaks/run_info.txt")
    # ls.append(output_path + "/peaks/all_treatments.igv.xml")
    return ls

rule bin_all:
    input:
        bin_targets