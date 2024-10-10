#' Normalize data to lie between 0 and 1
#'
#' Used for plotly_ridges
#' 
#' @param data the data frame to operate on
#' @param group optional column to group data on when normalizing
#' @param var_y the variable to become between 0 and 1
#' @param var_x the x variable
#' @param normalizer the function to do the normalizing
#' @param ycol when not grouping, an extra column to include
compute_group_normalization <- function(data,group,var_y,var_x,normalizer=function(x){(x-min(x))/(max(x)-min(x))},ycol=NULL){
    ycol = dplyr::select(data,{{ ycol }})
    ret = data |>
        dplyr::group_by( {{group}} ) |>
        dplyr::reframe(var_y=normalizer( {{var_y}} ),var_x= {{ var_x }},orig_y={{ var_y }} )
    cbind(ret,ycol)
}


#' Ridgeline Plot using Plotly
#'
#' Ridgeline plot using plotly. This allows combining different plots together and adjusting the scale via
#' a user specified normalizer function. 
#' 
#' @param data the dataframe to plot
#' @param x the x variable (e.g., length,age)
#' @param y the y variable (e.g., year)
#' @param height the ridgeline height variable (e.g., survey index)
#' @param xlab X axis label name (also used in hovertext)
#' @param ylab Y axis label name
#' @param heightlab name of the height variable (used in hovertext)
#' @param main main plot title
#' @param plot_group label for legend if plotting multiple
#' @param plot_color color of plot
#' @param normalizer the function to do the normalizing
#' @param n_group how to do the grouping for the normalizer
#' @param add prior plotly object if adding on to existing plot
#' @export
plotly_ridges <- function(data,x,y,height,xlab,ylab,heightlab,main,plot_group,plot_color,normalizer=function(x){(x-min(x))/(max(x)-min(x))},n_group=NULL,add=NULL){
    if(is.null(add)){
        pp = plotly::plot_ly()
    }else{
        pp = add
    }


    if(rlang::quo_is_null(enquo(n_group))){
                                        #normaled_dat = compute_group_normalization(data,group={{ n_group }},var_y={{ height }},var_x={{ x }},normalizer = normalizer)
        normaled_dat = compute_group_normalization(data,group=NULL,var_y={{ height }},var_x={{ x }},normalizer = normalizer,ycol={{ y }})
    }else{
        normaled_dat = compute_group_normalization(data,group={{ n_group }},var_y={{ height }},var_x={{ x }},normalizer = normalizer)
    }
    normaled_dat = normaled_dat[,c(unique(names(normaled_dat)))]

    fcolor = paste0("rgba(",paste0(col2rgb(plot_color),collapse = ","),",0.3)")

    syears = select(data,{{ y }} ) |>
        pull()
    syears = unique(syears)

    for(i in 1:length(syears)){
        ndat = normaled_dat |>
            dplyr::filter({{ y }} == syears[i])
        yy = ndat |>
            dplyr::mutate(yy = {{ y }} + var_y) |>
            dplyr::pull()

        xx = ndat |>
            dplyr::mutate(xx = var_x) |>
            dplyr::pull()

        og_y = ndat$orig_y
        df2 = data.frame(og_y=round(og_y,2),year=syears[i],og_x=xx)
        df2_list = split(df2,seq_len(nrow(df2)))

        hovtemp = paste0(xlab,": %{x}\n",heightlab,": %{customdata.og_y}\n",ylab,": %{customdata.year}")

        pp = pp |>
            plotly::add_trace(y=syears[i]+0.001,x=xx,mode="lines",line=list(color="white",width=0),showlegend=FALSE,hovertemplate="",hovermode=FALSE)
        if(i %in% c(1)){
            pp = pp |>
                plotly::add_trace(y=yy,x=xx,name=plot_group,lines=list(color=plot_color),mode="lines",fill="tonexty",legendgroup=plot_group,fillcolor=fcolor,showlegend=TRUE,customdata=df2_list,hovertemplate=hovtemp)
        }else{
            pp = pp |>
                plotly::add_trace(y=yy,x=xx,name=plot_group,lines=list(color=plot_color),mode="lines",fill="tonexty",legendgroup=plot_group,fillcolor=fcolor,showlegend=FALSE,customdata=df2_list,hovertemplate=hovtemp)
        }
    }


    pp = pp |>
        plotly::layout(yaxis=list(title=ylab),xaxis=list(title=xlab),title=main)
    pp
}



#' Normalizer function to put two different data sets on the same scale in a ridgeline plot.
#'
#' @param survey1 index column of first data set
#' @param survey2 index column of second data set
#' @param year optional column of years if you want to split data by years
#' @export
match_normalizer <- function(survey1,survey2,year=NULL){
    if(!is.null(year)){
        s1split = split(survey1,year)
        s2split = split(survey2,year)
        mins = unlist(mapply(function(x,y){
            rep(min(x,y),length(x))},x=s1split,y=s2split,SIMPLIFY=FALSE))
        maxs = unlist(mapply(function(x,y){
            rep(max(x,y),length(x))},x=s1split,y=s2split,SIMPLIFY = FALSE))

    }else{
        mins = rep(min(survey1,survey2),length(survey1))
        maxs = rep(max(survey1,survey2),length(survey1))
    }  
    cbind(mins=mins,maxs=maxs)
    mnorm <- function(x){
        (x-mins)/(maxs-mins)
    }
}

#' Normalizer function to put multiple different data sets on the same scale in a ridgeline plot
#'
#' @param ... index columns of the data sets
#' @param year optional column of years if you want to split data by years
#' @export
match_normalizer2 <- function(...,year=NULL){
    dots <- list(...)
    if(!is.null(year)){
        lsplits = lapply(dots,function(x){
            split(x,year)})
        mins = numeric(0)
        maxs = numeric(0)
        for(y in 1:length(unique(year))){
            cur_y = lapply(lsplits,function(x){x[[y]]})
            mins = c(mins,rep(do.call(min,cur_y),length(cur_y[[1]])))
            maxs = c(maxs,rep(do.call(max,cur_y),length(cur_y[[1]])))
        }
    }else{
        mins = rep(do.call(min,dots),length(...elt(1)))
        maxs = rep(do.call(max,dots),length(...elt(1)))
    }
    cbind(mins=mins,maxs=maxs)
    mnorm <- function(x){
        (x-mins)/(maxs-mins)
    }
}

    


#' Create a bubble plot using plotly
#'
#' @param data the data to plot
#' @param x the x variable column
#' @param y the y variable column
#' @param bub_z the variable to make into a bubble!
#' @param title the title of the plot
#' @param xlab the name of the x axis
#' @param ylab the name of the y axis
#' @export
plotly_bubbles <- function(data,x,y,bub_z,title,xlab,ylab){
    data = dplyr::mutate(data,symbol=19,
                  size = dplyr::case_when({{bub_z}} == 0 ~ 0.1,
                                   TRUE ~ abs({{bub_z}})),
                  clr = dplyr::case_when({{bub_z}} < 0 ~ "blue",
                                  TRUE ~ "red"),
                  value = dplyr::case_when({{bub_z}} < 0 ~ "negative",
                                    TRUE ~ "positive"))
   bubp = ggplot2::ggplot(data, aes(x={{x}}, y={{y}})) + 
	ggplot2::geom_point(aes(size = size,color=value),shape=data$symbol,alpha=0.6, show.legend = TRUE) +
       ggplot2::xlab(xlab) + ggplot2::ylab(ylab)+ggtitle(title) + ggplot2::scale_color_manual(values=c("blue","red"))
    p_bubp = plotly::ggplotly(bubp)
    p_bubp
}


#' plotly correlation matrix
#'
#' @param cov.matrix the covariance matrix to plot as correlation matrix
#' @param title Title of the plot
#' @param cov.names the names you want for cov.matrix (optional)
#' @export
plotly_corr_mat <- function(cov.matrix,title,cov.names=NULL){
    corr_mat = make_corr_mat(cov.matrix)
    if(!is.null(cov.names)){
        colnames(corr_mat) = cov.names
        rownames(corr_mat) = cov.names
    }

    cm = ggcorrplot::ggcorrplot(corr_mat) + ggplot2::ggtitle(title)
    pcm = plotly::ggplotly(cm) |>
        plotly::config(mathjax = "cdn")
    pcm
}

#' plotly residuals
#'
#' Plot residuals with horizontal line at zero and trendline
#'
#' @param data the data to be plotted
#' @param x the x variable
#' @param y the y variable
#' @param title the overall title
#' @param xlab the x axis label
#' @param ylab the y axis label
#' @export
plotly_resids <- function(data,x,y,title,xlab,ylab){
    r1 = ggplot2::ggplot(data) + ggplot2::geom_point(aes(x={{ x }},y={{ y }})) +
        ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::ggtitle(title) +
        ggplot2::geom_hline(yintercept=0) + ggplot2::stat_summary(aes(x={{ x }},y={{ y }}),fun=mean,geom="line",size=1,color="red")
    gr1 = plotly::ggplotly(r1)
    gr1
}


#' Basic Plot + ribbon
#'
#' @param x x variable to plot
#' @param y y variable to plot
#' @param ymin optional minimum y value of ribbon
#' @param ymax optional maximum y value of ribbon
#' @param plot_group the plot group of the plot to belong to
#' @param color the plot and ribbon color
#' @param title title of plot
#' @param xlab the x axis label of the plot
#' @param ylab the y axis label of the plot
#' @param add plot to add on to (the object)
#' @export
basic_plot <- function(x,y,ymin=NULL,ymax=NULL,
                       plot_group=NULL,color,title,xlab,ylab,add=NULL,CI_lab="95% CI"){
    if(is.null(add)){
        pp = plotly::plot_ly()
    }else{
        pp = add
    }

    fcolor = paste0("rgba(",paste0(col2rgb(color),collapse=","),"0.3)")
    f_pg = paste(plot_group,CI_lab)

    pp = pp |>
        plotly::add_trace(y=y,x=x,mode="lines",line=list(color=color),legendgroup=plot_group,name=plot_group)

    if(!is.null(ymin)){
        pp = pp |>
            plotly::add_ribbons(x=x,ymin=ymin,ymax=ymax,color=fcolor,legendgroup=f_pg,name=f_pg)
    }
    
    pp = pp |>
        plotly::layout(title=title,yaxis=list(title=ylab),xaxis=list(title=xlab))
    
        
    pp

}

#' Plot Probability Length at Age Matrix
#'
#' @param plamat the PLA at age matrix to plot
#' @param ages the vector of ages in the model
#' @param lengths the vector of lengths to use in the model
#' @param len_mus the vector of predicted mean lengths at age
#' @param len_sds the vector of predicted sd lengths at age
#' @param title the title for the plot
#' @param crit_value the critical value to use
#' @export
plotly_pla <- function(plamat,ages,lengths,len_mus,len_sds,title,crit_value=1.96){

    pladf = as.data.frame.table(plamat)
    levels(pladf$Var1) = ages
    levels(pladf$Var2) = lengths
    pladf$Var1 = as.numeric(as.character(pladf$Var1))
    pladf$Var2 = as.numeric(as.character(pladf$Var2))

    musdf = data.frame(age=ages,mu=len_mus) |>
        dplyr::mutate(upper = mu+crit_value*len_sds,
                      lower = mu-crit_value*len_sds)

    col = hcl.colors(12,"YlOrRd",rev=TRUE)
    
    plmp = ggplot2::ggplot(pladf,aes(Var1,Var2)) + ggplot2::geom_raster(aes(fill=Freq)) +
        ggplot2::xlab("Age") + ggplot2::ylab("Length Bin") +
        ggplot2::geom_line(aes(x=age,y=mu),data=musdf) +
        ggplot2::geom_line(aes(x=age,y=upper),musdf,linetype=2) +
        ggplot2::geom_line(aes(x=age,y=lower),musdf,linetype=2) +
        ggplot2::scale_fill_gradient(low=col[1],high=col[12]) +
        ggplot2::ggtitle(title)

    gplmp = plotly::ggplotly(plmp)
    gplmp
}
