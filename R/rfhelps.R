#' Soft max function
#'
#' @param x the vector to transform
#' @export
#'
soft_max_transform <- function(x){
    y = exp(x)/sum(exp(x))
    y
}

#' Soft max inverse
#'
#' Assumes c sums to zero
#'
#' @param y the vector to untransform
#' @export
inv_soft_max_transform <- function(y){
    z = log(y)
    x = z-mean(z)
    x
}




#' Transform monotone increasing vector to optimizer scale
#'
#' @param x the vector to transform
#' @export
ordered_transform <- function(x){
    y = numeric(length(x))
    y[1] = x[1]
    for(i in 2:length(x)){
        y[i] = log(x[i]-x[i-1])
    }
    y
}

#' Transform vector back to monotone increasing
#'
#' @param y the vector to transform back
#' @export
ordered_inv_transform <- function(y){
    x = numeric(length(y))
    x[1] = y[1]
    for(i in 2:length(y)){
        x[i] = x[i-1]+exp(y[i])
    }
    x
}

#' Logistic function with different maximum
#'
#' @param x the x parameter
#' @param L the maximum of the curve
#' @export
glogit <- function(x,L){
    L/(1+exp(-(x)))
}

#' Inverse logistic function with different maximum
#'
#' @param y the thing to inverse
#' @param L the maximum of the curve
#' @export
iglogit <- function(y,L){
    -log(L/y-1)
}

#' Simplex transform with any maximum
#'
#' @param x the vector to transform
#' @param L the maximum to hit
#' @export
simplex_transform <- function(x,L){
  zk =  x[-length(x)]
  K = length(x)

  for(i in 2:(length(zk))){
    zk[i] = x[i]/(L-sum(x[1:(i-1)]))
  }
  y = glogit(zk,L)-log(L/(K-1:(K-1)))
  y


}

#' Inverse simplex transform with any maximum
#'
#' @param the vector to inverse
#' @param L the maximum to hit
#' @export
inv_simplex_transform <- function(y,L){
  K = length(y)+1
  zk = iglogit(y+log(L/(K-1:(K-1))),L)
  x = zk
  for(i in 2:length(zk)){
    x[i] = (L-sum(x[1:(i-1)]))*zk[i]
  }
  x[K] = L-sum(x)
  x
}
    


#' Convert a vector that sums to 1 to be on all real space
#'
#' @param x the vector of porportions
#'
#' @export
unit_simplex_transform <- function(x){
  zk =  x[-length(x)]
  K = length(x)

  for(i in 2:(length(zk))){
    zk[i] = x[i]/(1-sum(x[1:(i-1)]))
  }
  y = qlogis(zk)-log(1/(K-1:(K-1)))
  y


}

#' Convert back to a vector that sums to 1
#'
#' @param y vector to convert back
#'
#' @export
inv_unit_simplex_transform <- function(y){
  K = length(y)+1
  zk = plogis(y+log(1/(K-1:(K-1))))
  x = zk
  for(i in 2:length(zk)){
    x[i] = (1-sum(x[1:(i-1)]))*zk[i]
  }
  x[K] = 1-sum(x)
  x
}


#' Create covariance matrix
#'
#' @param cov.fixed the covariance matrix (e.g from sd.report cov.fixed)
#' @param names optional vector of names to use
#' @export
make_corr_mat <- function(cov.fixed,names=NULL){
    xm = 0.5*(t(cov.fixed)+cov.fixed)
    if(!matrixcalc::is.positive.definite(xm)){
        warning("Not positive definite!")
    }
    xm1 = Matrix::nearPD(xm)$mat
    dsd = sqrt(Matrix::diag(xm1))
    corr.matrix = as.matrix(diag(1/dsd)%*%xm1%*%diag(1/dsd))
    if(is.null(names)){
        rownames(corr.matrix) = rownames(xm1)
        colnames(corr.matrix) = colnames(xm1)
    }else{
        rownames(corr.matrix) = names
        colnames(corr.matrix) = names
    }
    corr.matrix
}

#'Bin indices
#'
#' @param df the data.frame to bin
#' @param len.col the length column to bin
#' @param cols_to_bin columns to summarize
#' @param grouping extra columns to group by
#' @param breaks the breaks of the binning
#' @export
bin_indices <- function(df,len.col,cols_to_bin,grouping,breaks,...){
    df = df |>
        dplyr::mutate(binned=cut({{len.col}},breaks,include.lowest = TRUE)) |>
        dplyr::group_by(binned,{{grouping}}) |>
        dplyr::summarise(dplyr::across({{cols_to_bin}},sum))
    df
}


#' Extracts specified variable name estimate and std. error from an sdreport summary
#'
#' @param sdr The sdreport summary to pull from
#' @param var.name The variable to extract
#' @param ... Can specify dimensions as vectors, any number. MUST FOLLOW TEMPLATE ORDER & DIMENSIONS
#' @param flatten Convert array to flat data.frame
#' @export
extract_from_sdr <- function(sdr,var.name,...,flatten=FALSE){

    urnames = unique(row.names(sdr))
    if(!(var.name %in% urnames)){
        stop(paste(var.name,"not found in sdreport summary!"))
    }


    dats = sdr[row.names(sdr) == var.name,,drop=FALSE]
    est = dats[,1]
    std = dats[,2]


    ret = list(est=est,std=std)
    ret

    if(...length() > 0){
        dims = numeric(0)
        dimnameslist = list()
        dotnames = ...names()
        for(i in 1:...length()){
            dims = c(dims,length(...elt(i)))
            dimnameslist[[i]] = ...elt(i)
        }
        ret$est = array(ret$est,dim=dims,dimnames=dimnameslist)
        ret$std = array(ret$std,dim=dims,dimnames=dimnameslist)
        if(flatten){
            retE = as.data.frame.table(ret$est)
            temp = as.data.frame.table(ret$std)
            retE = cbind(retE,temp[,ncol(temp)])
            names(retE) = c(dotnames,paste0(var.name,".est"),paste0(var.name,".std"))
            ret = retE
        }

    }else{
        if(flatten){
            retE = data.frame(ret$est,ret$std)
            names(retE) = c(paste0(var.name,".est"),paste0(var.name,".std"))
            ret = retE
        }
    }


    ret
}

#' Create  confidence interval
#'
#' @param est thing of estimates
#' @param std thing of standard deviations
#' @param alpha value of alpha
#' @param undotrans function to undo transformation (e.g., exp)
#' @export
generate_CI <- function(est,std,alpha=0.05,undotrans=NULL){
    conf = 1-alpha/2
    crit = qnorm(conf)
    lower = est-crit*std
    upper = est+crit*std
    if(!is.null(undotrans)){
        lower = undotrans(lower)
        upper = undotrans(upper)
        est = undotrans(est)
    }

    list(lower=lower,upper=upper,est=est)
}

#'Convert to untransformed names
#'
#' Removes logs and ts and whatever from TMB names
#' @param var_names vector of variable names
#' @export
untransform_names <- function(var_names){
    tTMBN1 = gsub("^log_","",var_names)
    tTMBN2 = gsub("^t_","",tTMBN1)
    tTMBNtoo = gsub("engel\\.log_","engel.",tTMBN2)
    tTMBN3 = gsub("\\.log_","",tTMBNtoo)
    tTMBN4 = gsub("\\.t_","",tTMBN3)
    tTMBN4 = gsub("\\.logit_","",tTMBN4)
    tTMBN4 = gsub("^logit_","",tTMBN4)
    tTMBN4
}

#' Match TMB names with LaTeX or TeX names
#'
#' @param var_names vector of names to match
#' @param fixed boolean vector indicating whether or not a value is fixed
#' @param type either LaTeX or TeX
#' @param transformed use transformed form or not
#' @export
match_TMB_names <- function(var_names,fixed=rep(FALSE,length(var_names)),type,transformed){

    LaTeXN = c("\\(l\\)",
               "\\(L\\)",
               "\\(k\\)",
               "\\(s\\)",
               "\\(S\\)",
               "\\(\\sigma_{R}\\)",
               "\\(\\sigma_{N0}\\)",
               "\\(\\sigma_{Fa}\\)",
               "\\(\\sigma_{Fy}\\)",
               "\\(\\sigma_{I}\\)",
               "\\(q_{max,C}\\)",
               "\\(q_{50,C}\\)",
               "\\(q_{95,C}\\)",
               "\\(q_{max,E}\\)",
               "\\(q_{50,E}\\)",
               "\\(q_{95,E}\\)",
               "\\(\\rho_{EC}\\)",
               "\\(\\sigma_{T}\\)",
               "\\(\\sigma_{pe}\\)",
               "\\(\\phi_{M,y}\\)",
               "\\(\\phi_{M,a}\\)",
               "\\(\\sigma_{surv}\\)",
               "\\(\\sigma_{Q_max}\\)",
               "\\(K_{vb}\\)",
               "\\(L_{\\infty}\\)",
               "\\(t_{0}\\)",
               "\\(\\phi_{N,y}\\)",
               "\\(\\phi_{N,a}\\)",
               "\\(\\beta_{1,sel}\\)",
               "\\(\\beta_{2,sel}\\)",
               "\\(k_{sel}\\)",
               "\\(\\theta_{sel}\\)",
               "\\(\\sigma_{Cprop,a}\\)",
               "\\(\\rho_{Cprop}\\)",
               "\\(\\rho_{Q,C}\\)",
               "\\(\\bar{a}_{A,0}\\)",
               "\\(\\rho_I\\)"
               )

    lnames3 = gsub("\\\\\\(","",LaTeXN)
    lnames4 = gsub("\\\\\\)","",lnames3)
    TeXN = lapply(lnames4,plotly::TeX)

    ndf = data.frame(latex=LaTeXN,tex=unlist(TeXN))


    TMBN = c("log_ell",
             "log_L",
             "log_k",
             "log_s",
             "log_S",
             "log_recruit_sd",
             "log_N0_sd",
             "log_Fa_sd",
             "log_Fy_sd",
             "log_sd_survey",
             "log_Qmax",
             "log_QL50",
             "log_QL95",
             "engel.log_Qmax",
             "engel.log_QL50",
             "engel.log_QL95",
             "t_rho_ec",
             "log_landings_sd",
             "log_sd_pe",
             "t_phi_pe_year",
             "t_phi_pe_age",
             "log_surv_sd",
             "log_Q_max_sd",
             "growth.vbK",
             "growth.vbLinf",
             "growth.vbtzero",
             "t_phi_Npe_year",
             "t_phi_Npe_age",
             "log_b_beta1",
             "log_b_beta2",
             "log_sel_shape",
             "log_sel_scale",
             "log_sd_catch_prop_m",
             "logit_rhoC",
             "t_Q_rho_ar",
             "log_init_a_pg",
             "logit_rhoS")

    ## tTMBN1 = gsub("^log_","",TMBN)
    ## tTMBN2 = gsub("^t_","",tTMBN1)
    ## tTMBN3 = gsub("\\.log_","",tTMBN2)
    ## tTMBN4 = gsub("\\.t_","",tTMBN3)
    tTMBN4 = untransform_names(TMBN)
    

    if(transformed){
        rownames(ndf) = tTMBN4
    }else{
        rownames(ndf) = TMBN
    }
    if(type == "LaTeX"){
        ret = ndf[var_names,"latex"]
    }else{
        ret = ndf[var_names,"tex"]
    }
    ret = ifelse(fixed,paste0(ret,"*"),ret)
    ret
    
}               



#' Create Report
#'
#' Create a report from a model output file. Requires working Rmarkdown setup.
#'
#' @param title the title of the report and file name of the report
#' @param model_file the path to the rds output of the model
#' @param report_output where to output the report
#' @param tmb_data path of tmb input data
#' @param mod_data path of modified data
#' @param build build the html file?
#' @export
create_report <- function(title,model_file,report_output,tmb_data,mod_data,build=TRUE){
    
    pkgloc = system.file(package="rfhelpers")
    cmd = paste0(pkgloc,"/exec/report.Rmd")
###copy the report to the directory
    file.copy(cmd,to=paste0(report_output,"/","reportT.rmd"))
    parmys = list( modelpath=model_file,modelname=title,tmbdata=tmb_data,moddata=mod_data)
    outputdir = report_output
    outputfile = paste0(title,".html")
    rmarkdown::render(paste0(report_output,"/","reportT.rmd"),"html_document",output_file = outputfile,output_dir=outputdir,params=parmys)
    ##Remove the report
    file.remove(paste0(report_output,"/","reportT.rmd"))
}

#' Generate prior for conv. factor
#'
#' This uses the inverse of the conversion factor between
#' Engel and Campelen Trawl from Warren 1997
#' @param x the value of the length
#' @export
redfish_rat <- function(x){
    lny = 6.7580137+0.006839*x-1.927210*log(x)
    1/exp(lny)
}



#' Build model data
#'
#' Setup the model data
#'
#' @param years vector of years to include
#' @param ages vector of ages to use in the model
#' @param lengths vector of lengths to include
#' @param weight_array array of weights with dimnames for male/female, year, length, etc.
#' @param maturity_array array of maturity ogives to use, with appropriate dimnames
#' @param survey_df data.frame of survey data to use
#' @param landings_df data.frame of landings data, should include landings bounds
#' @param bin_ajdust midpoint of the length bins
#' @param min_size minimum required landing size
#' @param use_discards use discards or not
#' @param inf_length How far out to push the length matrix for the length plus group?
#' @param Q_prior_max Q prior max length
#' @export
build_data <- function(years,ages,lengths,weight_array,maturity_array,survey_df,landings_df,bin_adjust=0.5,min_size=22,use_discards=FALSE,inf_length=100,Q_prior_max=inf_length){
    ##In case of gaps
    Y = length(seq(years[1],years[length(years)]))
    A = length(ages)
    L = length(lengths)
    L3 =length(1:inf_length)

    start_length = lengths[1]
    end_length = lengths[length(lengths)]

    start_age = ages[1]
    end_age = ages[length(ages)]

    start_year = years[1]
    end_year = years[length(years)]

    weightsF = weight_array[as.character(1:inf_length),as.character(start_year:end_year),"female"]
    weightsM = weight_array[as.character(1:inf_length),as.character(start_year:end_year),"male"]

    maturityF = maturity_array[as.character(1:inf_length),as.character(start_year:end_year),"male"]
    maturityM = maturity_array[as.character(1:inf_length),as.character(start_year:end_year),"male"]

    survey2 = survey_df |>
        dplyr::filter(survey.year >= start_year) |>
        dplyr::filter(survey.year <= end_year) |>
        dplyr::mutate(length = dplyr::case_when(length <= start_length ~ start_length,
                                  length >= end_length ~ end_length,
                                  TRUE ~ length)) |>
        dplyr::group_by(survey.year,length) |>
        dplyr::mutate(total=sum(total),total.lcl=sum(total.lcl),total.ucl=sum(total.ucl),
               mean=sum(mean),mean.lcl=sum(mean.lcl),mean.ucl=sum(mean.ucl)) |>
        dplyr::distinct(survey.year,length,.keep_all=TRUE) |>
        dplyr::filter(total > 0)

    add.index <- function(data,scale=1){
        data$rsu = data$sample.units/mean(unique(data$sample.units))
        data$index = data$total/(scale*data$rsu)
        data
    }

    survey3 = add.index(survey2,1e6)
    
    survey_index = survey3$index
    survey_year = as.integer(survey3$survey.year-start_year)
    survey_length = as.integer(survey3$length-1)
    survey_type = ifelse(survey3$survey.year < 1995,0,1)

    landings_df = landings_df |>
        dplyr::filter(Years >= start_year) |>
        dplyr::filter(Years <= end_year)

    log_landingsL = log(landings_df$lower_bound)
    log_landingsU = log(landings_df$upper_bound)
    
    landing_year = as.integer(landings_df$Years-start_year)
    og_landings = landings_df$Total
    adj_landings = ifelse(is.na(landings_df$landadj),landings_df$Total,landings_df$landadj)
    discards = landings_df$discards

    if(use_discards){
        landing_nums = adj_landings
    }else{
        landing_nums = og_landings
    }

    min_size2 = rep(min_size,Y)

    agesindices = matrix(1:A,nrow=A,ncol=Y)
    vagesindices = as.vector(agesindices)
    yearsindices = matrix(1:Y,nrow=A,ncol=Y,byrow=TRUE)
    vyearsindices = as.vector(yearsindices)

    kages = matrix(agesindices,nrow=A*Y,ncol=A*Y,byrow=FALSE)
    kyears = matrix(yearsindices,nrow=A*Y,ncol=A*Y,byrow=FALSE)

    

    tmb.data = list(
        start_length = start_length,
        end_length = end_length,
        start_year = start_year,
        end_year = end_year,
        Y = Y,
        A = A,
        L = L,
        L2 = length(start_length:inf_length),
        L3 = L3,
        inf_length = inf_length,
        weightsF = weightsF,
        weightsM = weightsM,
        maturityF = maturityF,
        maturityM = maturityM,
        survey_index = survey_index,
        survey_year = survey_year,
        survey_length = survey_length,
        survey_type = survey_type,
        landing_year = landing_year,
        landing_nums = landing_nums,
        bin_adjust = bin_adjust,
        adj_landings = adj_landings,
        discards = discards,
        log_landingsL = log_landingsL,
        log_landingsU = log_landingsU,
        dum_indices = ifelse(is.na(discards),0,1),
        min_size = min_size2,
        Qpriors = c(redfish_rat(1:(Q_prior_max-1)),rep(redfish_rat(Q_prior_max),length(Q_prior_max:inf_length))),
        vagesindices = vagesindices-1,
        vyearsindices = vyearsindices-1,
        kages = kages-1,
        kyears = kyears-1
    )

    modDat = list(
        survey=survey3,
        landings=landings_df,
        discards = landings_df,
        weightsF=weightsF,
        weightsM=weightsM,
        maturityF=maturityF,
        maturityM=maturityM
    )
    
    ret = list(tmb.data=tmb.data,modified.data=modDat)
    ret

}





#'Run Model
#'
#' Run the model with the selected settings
#' @param weight_array array of weights with dimnames for male/female, year, length, etc.
#' @param maturity_array array of maturity ogives to use, with appropriate dimnames
#' @param survey_df data.frame of survey data to use
#' @param landings_df data.frame of landings data
#' @param F_type either 'age' or 'simple_blocked_len'
#' @param N_type either 'fixed' or 'cohort'
#' @param M_type Either 'fixed', 'M-devs', 'Lorenzen' or 'LorenzenDevs'
#' @param catch_type Either 'no_prop' or 'prop'
#' @param landing_type Either 'unbounded' or 'bounded'
#' @param base_M the base_M underlying value assumption
#' @param output_path path to put the report and model output/data (optional)
#' @param report_title what to call the report (optional)
#' @param upper_landing_bounds column of landings upper multiplier for bounds (only for bounded landings)
#' @param lower_landing_bounds column of landings lower multiplier for bounds (only for bounded landings)
#' @param catch_prop optional matrix with named row/cols of catch proportions
#' @param agg_key required if using catch_prop, specifies how to aggregate catch
#' @param years vector of years to include
#' @param ages vector of ages to use in the model
#' @param lengths vector of lengths to include
#' @param tmb.map optional named list, use to manually provide TMB map
#' @param random optional, override the default random effect parameters
#' @param control optional, override nlminb optimizer settings (see nlminb help)
#' @param start.parms optional named list of start parameters, used to override the starting parameters
#' @param data optional named list of data to give the TMB model, used to directly override the data
#' @param model_num optional model number for identity purposes
#' @param catch_sig_type Either 'independent', 'AR1' or 'compound'
#' @param survey_sig_type Either 'independent' or 'AR1'
#' @param inf_length the max length used in the model
#' @param optimize optimize the model?
#' @param use_sdreport run the sdreport?
#' @param plot_output generate the plot html?
#' @param use_Qpriors use the data based prior on survey Q
#' @param fix_L_inf if given a value in cm this will be used as the value of L_inf, setting L to the corresponding value
#' @param rand_N run the random N draw?
#' @param Q_prior_max the maximum of the Q prior length to use
#' @param inter_hessian use internal hessian instead of optimized
#' @param pg_ext controls how fine the plus group growth extension goes
#' @param gf_type controls some settings on the growth function 'noavgage' has no plus group extension, 'avgage' does.
#' @param inner.control inner.control to pass to MakeADFun
#' @param silent trace the model or not
#' @param rounding_bit how much 95% selectivity for catch is ahead of 50%
#' @param double_opt run the optmizer again for a bit?
#' @param dopt_control control for the second round of optimiztion
#' @param survey_sd_map list by year and length specifying how to internally map the survey sds
#' @param catch_prop_map list by year and length specifying how to internally map the catch prop sds
#' @export
#'
#' @useDynLib redfish

run_redfish_acl <- function(weight_array,maturity_array,survey_df,landings_df,F_type,N_type,M_type,catch_type,landing_type,base_M,
                            output_path=NULL,report_title=NULL,
                            catch_prop = NULL,
                            agg_key = NULL,upper_landing_bounds=NULL,
                            lower_landing_bounds=NULL, years=1983:2021,ages=1:20,lengths=7:45,tmb.map=NULL,random=NULL,control=NULL,start.parms=NULL,
                            data=NULL,
                            model_num=NULL,catch_sig_type="independent",survey_sig_type="independent",inf_length=60,optimize=TRUE,use_sdreport=TRUE,plot_output=FALSE,use_Qpriors=FALSE,fix_L_inf = NULL,rand_N=TRUE,Q_prior_max=35,inter_hessian=FALSE,pg_ext=200,gf_type="avgage",inner.control=list(smartsearch=TRUE),silent=FALSE,rounding_bit=0.075,double_opt=FALSE,dopt_control=NULL,survey_sd_map = NULL,catch_prop_map=NULL){

    orig_data = list(weight_array=weight_array,
                     maturity_array=maturity_array,
                     survey_df=survey_df,
                     landings_df=landings_df,
                     F_type=F_type,
                     N_type=N_type,
                     M_type=M_type,
                     catch_type=catch_type,
                     landing_type=landing_type,
                     base_M=base_M,
                     output_path = output_path,
                     report_title = report_title,
                     catch_prop = catch_prop,
                     agg_key = agg_key,
                     upper_landing_bounds = upper_landing_bounds,
                     lower_landing_bounds = lower_landing_bounds,
                     years=years,
                     ages = ages,
                     lengths = lengths,
                     tmb.map = tmb.map,
                     random=random,
                     control=control,
                     start.parms=start.parms,
                     data=data,
                     model_num=model_num,
                     catch_sig_type = catch_sig_type,
                     survey_sig_type = survey_sig_type,
                     inf_length = inf_length,
                     optimize = optimize,
                     use_sdreport=use_sdreport,
                     plot_output = plot_output,
                     use_Qpriors= use_Qpriors,
                     fix_L_inf = fix_L_inf,
                     rand_N = rand_N,
                     Q_prior_max=Q_prior_max,
                     inter_hessian = inter_hessian,
                     pg_ext = pg_ext,
                     gf_type = gf_type,
                     inner.control = inner.control,
                     silent = silent,
                     rounding_bit = rounding_bit,
                     double_opt = double_opt,
                     dopt_control = dopt_control,
                     survey_sd_map = survey_sd_map,
                     catch_prop_map = catch_prop_map
                     )
                     

    if(!is.null(lower_landing_bounds)){
        landings_df$lower_bound = lower_landing_bounds
        landings_df$upper_bound = upper_landing_bounds
    }else{
        landings_df$lower_bound = rep(0,nrow(landings_df))
        landings_df$upper_bound = rep(0,nrow(landings_df))
    }

        
    datas = build_data(years,ages,lengths,weight_array,maturity_array,survey_df,landings_df,inf_length=inf_length,Q_prior_max=Q_prior_max)

    tmb.data = datas$tmb.data

    if(N_type == "fixed"){
        tmb.data$N_type = 0
    }else{
        tmb.data$N_type = 1
    }

    if(F_type == "age"){
        tmb.data$F_type = 0
    }else if(F_type == "simple_blocked_len"){
        tmb.data$F_type = 5
    }
    
    
    
    
    
    
            
    tmb.data$GF_type = 3
    
    if(M_type == "fixed"){
        tmb.data$M_type = 0
    }else if(M_type == "M-devs"){
        tmb.data$M_type = 1
    }else if(M_type == "Lorenzen"){
        tmb.data$M_type = 3
    }else{
        tmb.data$M_type = 4
    }

    if(landing_type == "unbounded"){
        tmb.data$landing_type = 0
    }else{
        tmb.data$landing_type = 1
    }


    catch_list = list()
    if(!is.null(agg_key) & !is.null(catch_prop)){
        for(i in 1:ncol(agg_key)){
            temp_key = list()
            temp_key$year = as.numeric(colnames(agg_key)[i])-tmb.data$start_year
            temp_key$key = seq(1:tmb.data$L3)
            temp_key$key[temp_key$key < agg_key[1,i]] = agg_key[1,i]
            temp_key$key[temp_key$key > agg_key[2,i]] = agg_key[2,i]
            temp_key$key = temp_key$key-min(temp_key$key)
            temp_key$ysize = length(unique(temp_key$key))
            if(!is.null(catch_prop_map)){
                temp_key$mmap = catch_prop_map[[i]]
            }else{
                temp_key$mmap = rep(0,temp_key$ysize)
                names(temp_key$mmap) = agg_key[1,i]:agg_key[2,i]
            }
            propt = numeric(length(temp_key$key))
            propt[as.numeric(rownames(catch_prop))] = catch_prop[,colnames(agg_key)[i],drop=TRUE]
            propt2 = numeric(inf_length)
            propt2[1:inf_length] = propt[1:inf_length]
            propt2[inf_length] = sum(propt[inf_length:length(propt)])
            temp_key$prop = propt2
            catch_list[[i]] = temp_key
            names(catch_list)[i] = colnames(agg_key)[i]
        }
        ##Check if years in years
        CLyears = colnames(agg_key)[colnames(agg_key) %in% years]
        catch_list = catch_list[CLyears]
        
    }
    tmb.data$catch_list = catch_list

    tmb.data$rec_type = 0
    

    mapp = list()
    parms = list()
    randos = c("log_Fa","log_Fy")
 

    tmb.data$base_M = base_M
    ##Options not really used and kinda broken
    tmb.data$weight_type = 0
    tmb.data$cfg = list(sparse=FALSE,trace=TRUE)
    if(use_Qpriors == FALSE){
        tmb.data$Q_prior = 0
    }else{
        tmb.data$Q_prior = 1
        parms$log_Q_max_sd = log(0.5)
        ##parms$log_Q_rho = ordered_transform(log(seq(0.01,1,length.out = inf_len)))
        ##parms$log_Q_rho = rep(log(1),inf_len)
        ##parms$t_Q_rho_ar = rep(0.5,QAR_len)
        ##randos = c(randos,"log_Q_rho")
        mapp$log_Q_max_sd = as.factor(NA)
        ##mapp$log_Q_rho = as.factor(c(1:(Q_prior_max-1),rep(Q_prior_max,length(Q_prior_max:inf_length))))
        ##mapp$log_Q_rho = as.factor(c(1:tail(lengths,1),rep(tail(lengths,1),length((tail(lengths,1)+1):inf_length))))
    }
    
    
    tmb.data$Q_max_prior = 0.99

    ell_L = c(7,45)
    tran_ell_L = ordered_transform(log(ell_L))

    s_S = c(1,2)
    tran_s_S = ordered_transform(log(s_S))

    
    parms$log_ell = tran_ell_L[1]
    parms$log_L = tran_ell_L[2]
    parms$log_k = log(0.5)
    parms$log_S = s_S[2]
    parms$log_s = s_S[1]
    ##parms$log_S = tran_s_S[2]
    ##parms$log_s = tran_s_S[1]
    parms$log_recruit_sd = log(0.3)
    
    
    parms$log_N0_sd = log(0.3)
    parms$log_Fy = rep(log(0.2),tmb.data$Y)
    parms$log_Fy_sd = log(0.1)
    ##parms$log_sd_survey = log(0.2)
    ##parms$log_sd_survey_plus = log(0.2)
    ##parms$log_surv_plus_sd = log(0.28)
    parms$log_Qmax = log(1)
    parms$log_QL50 = log(12)
    parms$log_QL95 = log(20)
 
    parms$log_landings_sd = log(0.2)

       mapp$log_QL95 = as.factor(NA)
       mapp$log_QL50 = as.factor(NA)
       mapp$log_N0_sd = as.factor(NA)


    
    mapp$log_Fa = as.factor(c(NA,NA,1,2,3,rep(4,tmb.data$A-5)))
    mapp$log_landings_sd = as.factor(NA)
    mapp$log_Fa_sd = as.factor(NA)

    
    parms$engel.log_Qmax = log(1)
    parms$engel.log_QL50 = log(12)
    parms$engel.log_QL95 = log(20)



    if(M_type == "M-devs" | M_type == "LorenzenDevs"){
        parms$process_e = matrix(log(base_M),nrow=tmb.data$A,ncol=tmb.data$Y)
        parms$t_phi_pe_year = plogis(0.75)
        parms$t_phi_pe_age = plogis(0.75)
        parms$log_sd_pe = log(0.5)
        mapp$log_Qmax = as.factor(NA)
        randos = c(randos,"process_e")
    }

    if(N_type == "fixed"){
        ##parms$log_recruit = seq(5,1,length.out=tmb.data$Y)
        ##parms$log_recruit = rep(5,tmb.data$Y)
        parms$log_recruit = seq(5,5.1,length.out = tmb.data$Y)
        parms$log_N0 = rep(0.5,tmb.data$A-1)
        randos = c(randos,"log_recruit","log_N0")
    }else{
        parms$log_N_a = matrix(log(5),nrow=tmb.data$A,ncol=tmb.data$Y)
        parms$log_N_a[tmb.data$A,] = log(4.9)
        parms$log_N_a[1,] = seq(log(5),log(5.1),length.out = tmb.data$Y)
        parms$log_surv_sd = log(0.28)
        mapp$log_surv_sd = as.factor(NA)
        randos = c(randos,"log_N_a")
    }

    if(F_type == "age"){
        parms$log_Fa = c(-10,-10,rep(log(0.2),tmb.data$A-2))
        parms$log_Fa_sd = log(0.3)
    }else if(F_type == "simple_blocked_len"){
        tmb.data$mora_year = 1997-min(years)

        parms$log_b_beta1 = log(25)
        tmb.data$rounding_bit = rounding_bit
        ##parms$log_b_beta2 = log(35)

        parms$log_sel_scale = log(1)
        parms$log_sel_shape = log(15)
        randos = randos[randos != "log_Fa"]
        mapp$log_Fa = NULL
        mapp$log_Fa_sd = NULL
    }else if(F_type == "AR"){
        randos = randos[randos != "log_Fa"]
        randos = randos[randos != "log_Fy"]
        mapp$log_Fa = NULL
        mapp$log_Fa_sd = NULL

        parms$t_phi_F_age = qlogis(0.51)
        parms$t_phi_F_year = qlogis(0.51)
        parms$log_sd_F = qlogis(0.3)
        parms$array_F = matrix(log(0.3),nrow=tmb.data$A,ncol=tmb.data$L3)
        tmb.data$F_type = 3
        randos = c(randos,"array_F")
    }else if(F_type == "rwrw"){
        randos = randos[randos != "log_Fa"]
        mapp$log_Fa = NULL
        mapp$log_Fa_sd = NULL
        tmb.data$mora_year = 1997-min(years)
        

        parms$log_sd_rw1 = log(0.5)
        parms$log_sd_rw2 = log(0.5)
        parms$rw1 = rep(qlogis(0.3),tmb.data$L3)
        parms$rw2 = rep(qlogis(0.3),tmb.data$L3)
        randos = c(randos,"rw1","rw2")
        tmb.data$F_type = 4
    }
    
    
    
        


    if(catch_type == "prop"){
        if(is.null(catch_prop)){
            stop("catch_type prop needs supplied catch_prop")
        }else{
            tmb.data$catch_type = 1

            catch_map = lapply(catch_list,function(x){
                    x$mmap
            })
            catch_map = unlist(catch_map)
            
            parms$log_sd_catch_prop_m = rep(log(1),length(unique(catch_map)))
            ##parms$log_sd_catch_prop_m = matrix(log(1),nrow=length(catch_list[[1]]$key),ncol=length(catch_list))
            tkeys = do.call(cbind,lapply(catch_list,function(x){x$key}))
            ##mapp$log_sd_catch_prop_m = as.factor(tkeys)
            if(F_type == "age"){
                tmb.data$F_type = 0
            }
            
        }
        if(catch_sig_type == "independent"){
            tmb.data$catch_sig_sw = 0
        }else if(catch_sig_type == "AR1"){
            tmb.data$catch_sig_sw = 1
            parms$logit_rhoC = plogis(0.5)
        }else if(catch_sig_type == "compound"){
            tmb.data$catch_sig_sw = 2
            parms$logit_rhoC = plogis(0.5)
        }
        
        
    }else{
        tmb.data$catch_type = 0
    }

    if(!is.null(fix_L_inf)){
        tmb.data$GF_type = 2
        tmb.data$fixed_L_inf = fix_L_inf
        parms$log_L = NULL
    }
    

    if(landing_type == "bounded"){
        parms$log_landings_sd = log(0.02)
    }

    if(!is.null(start.parms)){
        parms[names(start.parms)] = start.parms
    }

    if(!is.null(tmb.map)){
        mapp[names(tmb.map)] = tmb.map
    }

    if(is.null(control)){
        control = list(iter.max=2000,eval.max=2000)
    }

    if(is.null(dopt_control)){
        dopt_control = control
    }
    

    if(!is.null(random)){
        randos = random
    }

    if(!is.null(data)){
        tmb.data[names(data)] = data
    }
    


    tmb.data$pg_ext = pg_ext

    if(gf_type == "avgage"){
        tmb.data$GF_type = 4
        parms$log_init_a_pg = log(pg_ext-pg_ext/2)
    }
    if(gf_type == 5){
        tmb.data$GF_type = 5
    }


    tmb.data$adr_flag = 0

    ##Creating the survey_list
    ## This makes life easier for AR survey indices

    s_ind = split(tmb.data$survey_index,tmb.data$survey_year)
    s_len = split(tmb.data$survey_length,tmb.data$survey_year)
    s_type = split(tmb.data$survey_type,tmb.data$survey_year)
    s_place = split(1:length(tmb.data$survey_index),tmb.data$survey_year)
    survey_list = list()
    for(i in 1:length(s_ind)){
        tlist = list()
        tlist$year = as.numeric(names(s_ind)[i])
        tlist$lengths = s_len[[i]]
        tlist$indices = s_ind[[i]]
        tlist$type = unique(s_type[[i]])
        tlist$place = s_place[[i]]-1
        tlist$nlen = length(tlist$indices)
        tlist$projy = 0
        if(is.null(survey_sd_map)){
            tlist$mmap = rep(0,length(s_ind[[i]]))
        }else{
            tlist$mmap = survey_sd_map[[i]]
        }
        
        survey_list[[i]] = tlist
    }
    tmb.data$survey_list = survey_list
    ## This is just the overall size of the survey indices because we need to recreate some vectors to not mess up plotting
    tmb.data$survey_size = length(tmb.data$survey_index)
    sur_map = lapply(survey_list,function(x){
        x$mmap
    })
    sur_map = unlist(sur_map)


    ##parms$log_sd_survey = log(0.5)
    parms$log_sd_survey = rep(log(0.5),length(unique(sur_map)))
    if(survey_sig_type == "AR1"){
        tmb.data$survey_sigtype = 1
        parms$logit_rhoS = 0.1
    }else{
        tmb.data$survey_sigtype = 0
    }

    ##
    tmb.data$r_proj = 0
    tmb.data$landing_proj_y = rep(0,length(tmb.data$landing_nums))
    tmb.data$cfg = list(sparse=FALSE,trace=TRUE)
    tmb.data$proj_type = 0
    tmb.data$proj_years = 0
    tmb.data$og_Y = tmb.data$Y
    tmb.data$supplied_F = rep(0.1,tmb.data$proj_years)
    
    obj = TMB::MakeADFun(tmb.data,parms,map=mapp,random=randos,intern=FALSE, DLL="redfish",inner.control = inner.control,silent = silent)
    
    
    if(optimize==TRUE){
        opt = nlminb(obj$par,obj$fn,obj$gr,control=control)
        if(double_opt == TRUE){
            opt = nlminb(opt$par,obj$fn,obj$gr,control=dopt_control)
            obj$fn()
        }
        

        if(inter_hessian == TRUE){
            finp = obj$env$parList(obj$env$last.par.best)
            jbo = TMB::MakeADFun(tmb.data,finp,map=mapp,random=randos,intern=TRUE, DLL="redfish",inner.control = inner.control,silent = silent)
            print("Calculating Internal Hessian...")
            hess = jbo$he()
            opt$internal_hess = hess
           }
 
    }else{
        opt = NA

        if(inter_hessian == TRUE){
            obj = TMB::MakeADFun(tmb.data,parms,map=mapp,random=randos,intern=TRUE, DLL="redfish",inner.control = inner.control,silent = silent)
        }
    }

    if(inter_hessian == FALSE){
        repp = obj$report()
    }else{
        if(optimize == TRUE){
            repp = obj$report()
        }else{
            repp = NA
        }
    }
    
    if(use_sdreport){
        obj$env$data$adr_flag = 1
        print("Performing sdreport")
        if(inter_hessian == FALSE){
            hess = optimHess(opt$par,obj$fn,obj$gr)
            sdr = TMB::sdreport(obj,hessian.fixed = hess)
        }else{
            sdr = TMB::sdreport(obj,hessian.fixed = hess)
        }
        
        ssdr = summary(sdr)

        if(N_type == "cohort"){
            indN = which(names(sdr$value) == "resN")
            siggy = sdr$cov[indN,indN]
            tryCatch({
            if(rand_N){
                resND = MASS::mvrnorm(1,sdr$value[indN],sdr$cov[indN,indN])
                repp$resND = resND
                    
                
            }})
        }else{
            siggy = NA
        }

        ## if(catch_type == "prop"){
        ##     indC = grep("std_diffs_*",names(sdr$value))
        ##     resCD = MASS::mvrnorm(1,sdr$value[indC],sdr$cov[indC,indC])
        ##     repp$resCD = resCD
        ## }
        
        fixed_names = unique(names(opt$par))
        untrans_names = untransform_names(fixed_names)

        cov.trans = sdr$cov[names(sdr$value) %in% untrans_names,names(sdr$value) %in% untrans_names]
        rownames(cov.trans) = names(sdr$value)[names(sdr$value) %in% untrans_names]
        sdr$cov <- NULL
    }else{
        sdr = NA
        ssdr = NA
        cov.trans = NA
        siggy = NA
        hess = NA
    }
    

    keeps = list(obj=obj,map=mapp,sparms=parms,randos=randos,opt=opt,report=repp,sdreport=sdr,ssdr=ssdr,cov.trans=cov.trans,call=match.call(),orig_data=orig_data,sigN=siggy,hessian=hess)

    
    if(!is.null(output_path)){
        if(is.null(model_num)){
            model_num = ""
        }

        tmbdpathS = paste0(output_path,report_title,model_num,"tmbdata.rds")
        moddpathS = paste0(output_path,report_title,model_num,"moddata.rds")
        keepspathS = paste0(output_path,report_title,model_num,"out.rds")

        if(!dir.exists(output_path)){
            dir.create(output_path)
        }
        
        
        saveRDS(tmb.data,tmbdpathS)
        saveRDS(datas$modified.data,moddpathS)
        saveRDS(keeps,keepspathS)
        if(plot_output == TRUE){
            ##Rmarkdown works relative to the path of where the report is going
            tmbdpath = paste0("./",report_title,model_num,"tmbdata.rds")
            moddpath = paste0("./",report_title,model_num,"moddata.rds")
            keepspath = paste0("./",report_title,model_num,"out.rds")

            create_report(report_title,
                          model_file=keepspath,
                          report_output=output_path,
                          tmb_data=tmbdpath,
                          mod_data=moddpath)
        }
    }

    ret = list(output=keeps,tmb.data=tmb.data,mod.data=datas$modified.data,orig_data=orig_data)
    if(optimize == TRUE){
        print(keeps$opt$convergence)
    }
    ret
    
    
}

#' Run Projections on the model
#'
#' Project the model based on certain assumptions. If both supplied_F and given_catch are null we
#' assume to just use the terminal year. Projections are currently only implemented for models using
#' 'simple_blocked_len' selectivity and those with cohort effects. 
#' 
#'
#' @param model.object the output of run_redfish_acl
#' @param proj_years number of years to project
#' @param weights the weights at length for the projection years
#' @param maturity the maturity ogives at length for the projection years
#' @param supplied_F optional supplied vector to use for log F_y
#' @param given_catch optional vector of catch to use for projections
#' @export
run_projections <- function(model.object,proj_years=3,weights,maturity,supplied_F=NULL,given_catch=NULL){
    tdat = model.object$tmb.data
    gparms = as.list(model.object$output$sdreport,"Estimate")
    ttmap = model.object$output$map
    tdat$adr_flag = 1
    
    

    if(is.null(supplied_F) & is.null(given_catch)){
 
    
        if(!is.null(gparms$log_N_a)){
            oldNalen = length(gparms$log_N_a)
            lastCol = gparms$log_N_a[,ncol(gparms$log_N_a)]
            gparms$log_N_a =  cbind(gparms$log_N_a,matrix(rep(lastCol,proj_years),nrow=model.object$tmb.data$A,proj_years))
        }
    
        supplied_F = rep(gparms$log_Fy[model.object$tmb.data$Y],proj_years)
        goose = as.character(ttmap$log_Fy)
        goose[(model.object$tmb.data$Y+1):((model.object$tmb.data$Y+proj_years))] = (model.object$tmb.data$Y+1):((model.object$tmb.data$Y+proj_years))
        goose = as.factor(goose)
        ttmap$log_Fy = goose
    
    

        gparms$log_Fy = c(gparms$log_Fy,supplied_F)
        tdat$orig_end_year = tdat$end_year
        tdat$end_year = tdat$end_year+proj_years
        tdat$Y = tdat$Y+proj_years
        ##Yes this is on purpose
        tdat$og_Y = tdat$Y
        tdat$maturityF = cbind(tdat$maturityF,maturity)
        tdat$weightsF = cbind(tdat$weightsF,weights)
        tdat$r_proj = 1
    
        obj2 = TMB::MakeADFun(tdat,gparms,random=model.object$output$randos,map=ttmap,DLL="redfish")
        obj2$fn()
        rep2 = obj2$report()
        sdr2 = TMB::sdreport(obj2,hessian.fixed = model.object$output$hessian)

        keeps = list(obj=obj2,map=ttmap,sparms=gparms,randos=model.object$output$randos,opt=model.object$output$opt,report=rep2,sdreport=sdr2,ssdr=summary(sdr2),cov.trans=model.object$output$cov.trans,call=match.call(),orig_data=model.object$output$orig_data,sigN=model.object$output$sigN,hessian=model.object$output$hessian)

        ret = list(output=keeps,tmb.data=tdat,mod.data=model.object$mod.data,orig_data=model.object$orig_data)

    
    }else{

          if(!is.null(gparms$log_N_a)){
            oldNalen = length(gparms$log_N_a)
            ##gparms$log_N_a =  cbind(gparms$log_N_a,matrix(0,nrow=model.object$tmb.data$A,proj_years))
            lastCol = gparms$log_N_a[,ncol(gparms$log_N_a)]
            gparms$log_N_a =  cbind(gparms$log_N_a,matrix(rep(lastCol,proj_years),nrow=model.object$tmb.data$A,proj_years))
        }
    
    
    

        tdat$orig_end_year = tdat$end_year
        tdat$end_year = tdat$end_year+proj_years
        tdat$og_Y = tdat$Y
        ##Again on purpose
        tdat$Y = tdat$Y+proj_years
        tdat$maturityF = cbind(tdat$maturityF,maturity)
        tdat$weightsF = cbind(tdat$weightsF,weights)
        tdat$r_proj = 1
        if(!is.null(supplied_F)){
            tdat$proj_type = 0
            tdat$supplied_F = supplied_F
            randos = model.object$output$randos
        }else{
            tdat$proj_type = 1
            tdat$given_catch = given_catch
            ##gparms$proj_Fys = rep(log(0.1),proj_years)
            randos = model.object$output$randos
            ##randos = c(model.object$output$randos,"proj_Fys")
        }
        
    
        obj2 = TMB::MakeADFun(tdat,gparms,random=randos,map=ttmap,DLL="redfish")
        obj2$fn()
        rep2 = obj2$report()
        sdr2 = TMB::sdreport(obj2,hessian.fixed = model.object$output$hessian)

        keeps = list(obj=obj2,map=ttmap,sparms=gparms,randos=model.object$output$randos,opt=model.object$output$opt,report=rep2,sdreport=sdr2,ssdr=summary(sdr2),cov.trans=model.object$output$cov.trans,call=match.call(),orig_data=model.object$output$orig_data,sigN=model.object$output$sigN,hessian=model.object$output$hessian)

        ret = list(output=keeps,tmb.data=tdat,mod.data=model.object$mod.data,orig_data=model.object$orig_data)


        
    }

    ##Check if random effects have varied greatly from the model fit and if so throw a warning
    proFy = ret$output$report$log_Fy[1:model.object$tmb.data$Y]
    modFy = model.object$output$report$log_Fy
    FyC = all.equal(proFy,modFy)
    if(FyC != TRUE){
        warning("F Random Effects have changed from model fit!")
    }
    
    proN = ret$output$report$N[,1:model.object$tmb.data$Y]
    modN = model.object$output$report$N
    NyC = all.equal(proN,modN)
    if(NyC != TRUE){
        warning("N Random Effects have changed from model fit!")
    }
    
    

    
    ret
}


#' Run retrospective peels
#'
#' Construct models removing the number of specified years. Returns a list of model objects for each of the last years removed.
#' 
#' 
#' @param model.object the output of run_redfish_acl
#' @param n_peels the number of years to peel back, including the original
#' @param use_sdr use sdreport or not?
#' @param quick_start use fixed effect parameters from original model fit as starting values or use defaults (false)
#' @param silence disable TMB tracing for each peel?
#' @export
run_peels <- function(model.object,n_peels=6,use_sdr=TRUE,quick_start=TRUE,silence = TRUE){
    od = model.object$orig_data

    qpars = as.list(model.object$output$sdr,"Estimate")
    qpars = qpars[!(names(qpars) %in% model.object$output$randos)]
    if(quick_start == FALSE){
        qpars = od$start_parms
    }

    peels = list()
    peels[[1]] = model.object
    for(i in 2:n_peels){
        print(paste("On peel:",i))
        yearsP = head(od$years,length(od$years)-(i-1))

        peels[[i]] = with(od,run_redfish_acl(weight_array = weight_array,maturity_array = maturity_array,
                                             survey_df = survey_df,landings_df = landings_df,N_type = N_type,
                                             M_type = M_type,F_type=F_type,landing_type = landing_type,catch_type = catch_type,
                                             base_M = base_M,output_path = NULL,report_title = NULL,catch_prop=catch_prop,agg_key=agg_key,
                                             upper_landing_bounds = upper_landing_bounds,
                                             lower_landing_bounds = lower_landing_bounds,
                                             years=yearsP,ages= ages,
                                             lengths=lengths,tmb.map=tmb.map,
                                             random=random,control=control,start.parms = qpars,
                                             data=data,model_num=model_num,catch_sig_type = catch_sig_type,
                                             survey_sig_type = survey_sig_type,
                                             inf_length = inf_length,
                                             optimize=TRUE,use_sdreport = use_sdr,
                                             plot_output=FALSE,use_Qpriors=use_Qpriors,fix_L_inf=fix_L_inf,
                                             rand_N=rand_N,Q_prior_max=Q_prior_max,inter_hessian=inter_hessian,
                                             pg_ext=pg_ext,gf_type=gf_type,inner.control=inner.control,
                                             silent=silence,rounding_bit=rounding_bit,double_opt=double_opt,
                                             survey_sd_map = survey_sd_map,catch_prop_map))

        ##Update starting pars to last peel for hopefully faster finishes?
        qparsT = as.list(peels[[i]]$output$sdr,"Estimate")
        qparsT = qparsT[!(names(qparsT) %in% model.object$output$randos)]
        if(quick_start == TRUE){
            qpars = qparsT
        }
    }
    peels
}

                          
#' Calculate Mohn's rho and matrix for a specified quantity
#'
#' 
#' @param peels.object the object containing all the peels
#' @param quantity name of quantity to get
#' @param details report the details?
#' @param get_years grab the years to use?
#' @param ... things to pass to extract_from_sdr, if using years call it year
#'
#' @export
get_mohns <- function(peels.object,quantity,details=FALSE,get_years=TRUE,...){

    quant_list = list()
    dottys = list(...)

    for(i in 1:length(peels.object)){
        ssdr = peels.object[[i]]$output$ssdr
        if(get_years == FALSE){
            arg_list = list(sdr = ssdr, var.name = quantity, flatten = TRUE)
        }else{
            tyear = peels.object[[i]]$orig_data$years
            arg_list = list(sdr = ssdr, var.name = quantity, flatten = TRUE,year=tyear)
        }
        arg_list = c(arg_list, dottys)
        if (!is.null(dottys$year)) {
            arg_list$year = head(dottys$year, length(dottys$year) - 
                i)
        }
        temp = do.call(extract_from_sdr, arg_list)
        names(temp)[grep("\\.est", names(temp))] = paste0(names(temp)[grep("\\.est", 
                                                                           names(temp))], "_peel", i)
        names(temp)[grep("\\.std", names(temp))] = paste0(names(temp)[grep("\\.std", 
                                                                           names(temp))], "_peel", i)
        quant_list[[i]] = temp
    }

    quant_combo = dplyr::left_join(quant_list[[1]], quant_list[[2]])
    for (i in 3:length(quant_list)) {
        quant_combo = dplyr::left_join(quant_combo, quant_list[[i]])
    }
    estimates = dplyr::select(quant_combo, dplyr::contains(".est"))
    mohns = icesAdvice::mohn(estimates, peels = length(peels.object) - 
        1, details = details)
    ret = list(mohn = mohns, data = quant_combo,estimate=estimates)
    ret

}



#' Convert from schnute VB k to traditional VB K
#' @param k the k parameter value of the Schnute VB growth curve
#' @export
#'
get_K <- function(k){
    K = -log(k)
    K
}

#' Get L_infty from Schnute VB
#' @param L the largest length at age mean
#' @param ell the smallest length at age mean
#' @param k the Schnute k parameter
#' @param M the number of age classes
#'
#' @export
get_L_infinity <- function(L,ell,k,M){
    top = L-ell*k^(M-1)
    bot = 1-k^(M-1)
    Linf = top/bot
    Linf
}

#' Get t zero from Schnute VB
#' @param L the largest length at age mean
#' @param ell the smallest length at age mean
#' @param k the Schnute k parameter
#' @param M the number of age classes
#' @param a1 first age class in model
#'
#' @export
get_t_zero <- function(L,ell,k,M,a1){
    part = (L-ell)/(L-ell*k^(M-1))
    t_zero = a1-(1/log(k))*log(part)
    t_zero
}


#' Retro Report
#'
#' Create a report for the retros run, plots SSB, Average F, survey Abundance, abundance at age
#'
#' @param peels.object the object containing the peels
#' @param report_output the path to put the report
#' @param title name of the file to make
#' @export
retro_report <- function(peels.object,report_output,title="RetroReport"){
    
    pkgloc = system.file(package="rfhelpers")
    cmd = paste0(pkgloc,"/exec/retros.Rmd")
###copy the report to the directory
    file.copy(cmd,to=paste0(report_output,"/","retrosT.rmd"))
    parmys = list(peelsobj=peels.object)
    outputdir = report_output
    outputfile = paste0(title,".html")
    rmarkdown::render(paste0(report_output,"/","retrosT.rmd"),"html_document",output_file = outputfile,output_dir=outputdir,params=parmys)
    ##Remove the report
    file.remove(paste0(report_output,"/","retrosT.rmd"))
}

#' Projection Report
#'
#' Create a report for a set of projections
#'
#' @param proj.object the object containing the projections
#' @param report_ouput the path to put the report
#' @param title name of the file to make
#' @export
proj_report <- function(proj.object,report_output,title="ProjReport"){
    pkgloc = system.file(package="rfhelpers")
    cmd = paste0(pkgloc,"/exec/projections.Rmd")
###copy the report to the directory
    file.copy(cmd,to=paste0(report_output,"/","projectionsT.rmd"))
    parmys = list(projobj=proj.object)
    outputdir = report_output
    outputfile = paste0(title,".html")
    rmarkdown::render(paste0(report_output,"/","projectionsT.rmd"),"html_document",output_file = outputfile,output_dir=outputdir,params=parmys)
    ##Remove the report
    file.remove(paste0(report_output,"/","projectionsT.rmd"))   
}
