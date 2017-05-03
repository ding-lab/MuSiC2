# MuSiC2 ClinicalCorrelation

# by Qunyuan Zhang (qunyuan@wustl.edu)
# Matthew Wyczalkowski (m.wyczalkowski@wustl.edu)

# Usage: 
#   Rscript ClinicalCorrelation.pm [-v][-d] model.file variant.file trait.file out.file 

# This program is for Generalized Linear Model (GLM) analysis
# Y = X + covariate_1 + covariate_2 ......
# 
# Y is a quantitative or categorical trait (e.g., clinical trait)
# X is a quantitative or categorical variant (e.g., gene or mutation)

# Canonical Names:
# X - variant.  Aka matrix_file
# Y - trait.  Aka clinical_data

options(warn=1) 
options("width"=200) # DEBUG


# Return the command line argument associated with a given flag (i.e., -o foo),
# or the default value if argument not specified.
# Note that this will break if an argument is not supplied after the flag.
get_val_arg = function(args, flag, default) {
    ix = pmatch(flag, args)
    if (!is.na(ix)){ val = args[ix+1] } else { val = default }
    return(val)
}

# Return boolean specifying whether given flag appears in command line (i.e., -o),
get_bool_arg = function(args, flag) {
    ix = pmatch(flag, args)
    if (!is.na(ix)){ val = TRUE } else { val = FALSE }
    return(val)
}

# Usage:
#   args = parse.args()
#   print(args$disease.filter)
parse.args = function() {
    args = commandArgs(trailingOnly = TRUE)

    # optional arguments
    verbose = get_bool_arg(args, "-v")
    debug = get_bool_arg(args, "-d")

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    out.file = args[length(args)];      args = args[-length(args)]
    trait.file = args[length(args)];    args = args[-length(args)]
    variant.file = args[length(args)];  args = args[-length(args)]
    model.file = args[length(args)];    args = args[-length(args)]

    val = list('verbose'=verbose, 'debug'=debug, 'out.file'=out.file, 
               'variant.file'=variant.file, 'trait.file'=trait.file, 'model.file'=model.file)
    if (val$verbose) { print(val) }

    return (val)
}


# model.file is TSV with columns:
# "type":  Q=quantitative trait, B=binary trait
# "model.trait.name":     trait name
# "model.variant.name.string":     variant/gene name; if x=NA or blank, it will be determined by variant.file and all.variant.names
# "model.covariate":  covariate(s)
# "memo":  arbitrary text 

# trait.file is TSV with column headers, column 1 sample ID

# variant.file is TSV data file 
# the first column must be sample id (the same as in trait.file)
# must be defined 

# Multiple variants can be deefined in model with, model.variant.name.string="var_1|var_2|var_3"
# if model.variant.name.string *, "NA", or blank, will use all column names in variant.file as x variable names
# these must be found in column names of variant.file  

# This will crash if trait is all NA
do.glm=function(merged.data,trait,variant,covariate=NA,analysis.type) {

    if (!is.na(covariate)) {
        if (length(levels(merged.data[,c(covariate)])) == 1) {
            warning("One level in covariate detected.  Not doing covariate.");
            covariate = NA
        }
    }
#    print(merged.data[,c(variant, trait)])
#    print(covariate)

    if (nchar(covariate)==0 | is.na(covariate) | is.null(covariate)) { 
        model=formula(paste(trait,"~",variant)) 
    } else {
        model=formula(paste(trait,"~",covariate,"+",variant))
    }

    if (analysis.type=="B") 
        family=binomial(link = "logit")
    if (analysis.type=="Q") 
        family=gaussian(link = "identity")
    fit=glm(formula=model,data=merged.data,family=family)
    return(fit)
}

get.merged.variant.trait.data = function(trait.file, variant.file) {
    trait.data = read.table(trait.file,na.strings = c("","NA"),sep="\t",header=T)

    # TODO: Here, it would be smart to check whether all model names (model[,2]) correspond to columns in y.
    # obscure errors result if people mistype or use remapped characters (:-<=>)

    variant.data = read.table(variant.file,na.strings = c("","NA"),sep="\t",header=T)
    variant.id = colnames(variant.data)[1]
    all.variant.names = colnames(variant.data)[-1]  # all variant data column names except the first

    trait.id = colnames(trait.data)[1]
    ysid = !(colnames(trait.data) %in% all.variant.names)  # ysid are trait names not in variant names
    trait.data = trait.data[,ysid]                     # Subset trait data to those columns not in variant data (?)

    if (sum(ysid)==1) {  # ???  This seems to imply complete overlap in column names between variant and trait
        trait.data = data.frame(id=trait.data)
        colnames(trait.data)[1] = trait.id
    }

    # Variant and trait data merged here
    merged.data=merge(variant.data,trait.data,by.x = variant.id, by.y = trait.id)

    return( list('merged.data'=merged.data, 'all.variant.names'=all.variant.names) )
}

run.analysis = function(model, merged.data, all.variant.names, debug=FALSE) {
    glm.results=NULL

    for (i in c(1:nrow(model))) {  # loop through all rows in model  TODO: rewrite using apply()
    # analysis_type   clinical_data_trait_name    variant/gene_name   covariates  memo
        analysis.type=model[i,1]
        model.trait.name=model[i,2]
        model.variant.name.string=model[i,3]
        model.covariate=model[i,4]
        memo=model[i,5]

        model.variant.names = parse.variant.name.string(model.variant.name.string, all.variant.names)

        if (length(model.covariate)==0) 
            model.covariate=NA

        for (variant.name in model.variant.names) {
            if (debug) {
                cat(paste("    Processing: trait =", model.trait.name, " variant =", variant.name, "analysis =", analysis.type, " covariate =", model.covariate, "\n") )
            }
            result = evaluate.glm(analysis.type, merged.data, model.trait.name, variant.name, model.covariate) 
            if (!is.null(result)) {
                glm.results = rbind(glm.results, cbind(model.trait.name, analysis.type, variant.name, result$fit.result, result$coeff, model.covariate, memo))
            }
        }
    }

#"model.trait.name","analysis.type","variant.name","Df","Deviance","Resid. Df","Resid. Dev","F","Pr(>F)","model.covariate","memo"
    colnames(glm.results) = c("y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance", 
                              "F_statistic","p-value","coefficient","covariates","memo")
    glm.results$FDR = p.adjust(glm.results[,"p-value"], method="fdr") # Benjamini, Y., and Hochberg, Y. (1995). http://www.jstor.org/stable/2346101

    return(glm.results)
}

parse.variant.name.string = function(model.variant.name.string, all.variant.names) {
    if (!is.na(model.variant.name.string) & nchar(model.variant.name.string)>0) 
        model.variant.names=strsplit(model.variant.name.string,split="[|]")[[1]]
    if (is.na(model.variant.name.string)[1] | nchar(model.variant.name.string)[1]==0 | model.variant.name.string=="*") 
        model.variant.names=all.variant.names 
    return(model.variant.names)
}

evaluate.glm = function(analysis.type, merged.data, model.trait.name, variant.name, model.covariate) {
# Returns result of glm analysis.  Columns are standardized so that B and Q return row of same column names

    glm = try(do.glm(merged.data,model.trait.name,variant.name,model.covariate,analysis.type)) 

    # Type of error to catch here:
    #   Error in family$linkfun(mustart) : 
    #     Argument mu must be a nonempty numeric vector
    # This occurs if too many NA's
    # This also occurs: Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : contrasts can be applied only to factors with 2 or more levels
    #    above happens when using covariate and e.g. disease BLCA has all NA even if other diseases do not

    if(class(glm)[1] == "try-error") {
        cat(paste("    Error caught in GLM, continuing.  model.trait.name =", model.trait.name, 
                  " variant.name =", variant.name, " model.covariate =", model.covariate, "\n") )
        return(NULL);
    }

    if (analysis.type=="Q") {
        test = "F"
    } else if (analysis.type=="B") {
        test = "Chisq"
    } else {
        stop("Unknown model analysis.type ", analysis.type) 
    }

    fit = try(anova(glm,test=test))
    coeff = coefficients(glm)[[variant.name]] 
    if (class(fit)[1]!="try-error") {
        fit=as.matrix(fit)

        # Column names for Q: "Df"         "Deviance"   "Resid. Df"  "Resid. Dev" "F"          "Pr(>F)"
        # Column names for B: "Df"         "Deviance"   "Resid. Df"  "Resid. Dev"           "Pr(>Chi)"
        # To allow B and Q data to be mixed, we:
        # 1. Add an "F" column with value "NA" to any B result
        # 2. Rename last column to just "Pr"

        if (variant.name %in% rownames(fit)) {
            fit.result = as.data.frame(t(fit[variant.name,]))
            if (analysis.type == "B") {
                # want: Df Deviance "Resid. Df" "Resid. Dev" F "Pr(>F)"
                fit.result$F=NA
                fit.result=fit.result[,c(1:4,6,5)]
            }
            colnames(fit.result)[6] = "Pr"

            glm.result = list('fit.result'=fit.result, 'coeff'=coeff)

            # if this is B, add $F_statistic
            return(glm.result)
        }
    } else {
        cat(paste("    Error caught in anova, continuing.  model.trait.name =", model.trait.name, 
                  " variant.name =", variant.name, " model.covariate =", model.covariate, "\n") )
    }
    return(NULL);
}


args = parse.args()

model = read.table(args$model.file,colClasses="character",na.strings = c("","NA"),sep="\t",header=T)
data = get.merged.variant.trait.data(args$trait.file, args$variant.file)
glm.results = run.analysis(model, data$merged.data, data$all.variant.names, args$debug)

if (!is.null(glm.results)) {
    write.table(glm.results,file=args$out.file,quote=F,sep="\t",row.names=F) 
}
