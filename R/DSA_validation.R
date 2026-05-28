#' Use DSA algorithm to get optimal grouping result for the target country
#'
#' @param target target country for stepwise model
#' @param Pool Candidate Pool
#' @param COG is the maximum group size when to stop the algorith.
#' @param N maximum number of country-specific components
#' @param LossFunction SSE / MAE / MSPE / MAPE
#'
#' @return The best Pool, COG, N
#' @export
#'
DSA_V5<-function(target,Pool,COG,N,LossFunction)
  ## Pool: Candidate Pool
  ## PoolStatus: Record whether countries in the pool is available to be considered.(0-not available since it has been in the group; 1-available for current move)
  ## the criterion we use would be in sample SSE (each time choose one with smallest in sample SSE)
  ## V2: add N as tuning parameter
  ## V5: Treat female and male equally as populations
{
  Best=matrix(list(),nrow=N, ncol=COG) ## record optimal grouping for each size
  LOSS=matrix(0,nrow=N,ncol=COG) ## record in-sample MAPE of target country for all optimal groups and number of country-specific components using revised model
  Poolnames=names(Pool) ## record the original candidate pool
  #PoolStatus=as.vector(rep(1,length(Pool)))
  PoolStatus=as.vector(rep(1,length(Pool)))
  ## When group size=1
  ## the group is target itself
  ## the relevant SSE would be the in-sample SSE with component_size country-specific components
  for(component_size in 1:N)
  {
    Best[[component_size,1]]= list(Pool[[which(names(Pool)==target$label)]])#target
    LOSS[component_size,1]=loss(prefit_v2(target,target,component_size,Pool))[[LossFunction]]
  }
  for(component_size in 1:N)
    ## do D/S/A with different number of country-specific components
  {
    size=1 ## group size, when size=1, only the target country in the group
    PoolStatus=as.vector(rep(1,length(Pool)))
    ## Start with the first addition
    error=Inf
    i=1
    for(country in Pool)
    {
      label=country$label
      temp=Groupdemogdata_V2(Pool[[which(names(Pool)==target$label)]],########### used to be target
                          country,country$year[1],country$year[length(country$year)])
      error_new=loss(prefit_v2(temp,Target=target,component_size,Pool))[[LossFunction]]
      if(error_new<error)
      {
        error=error_new
        Choice=label
        Index=i
      }
      i=i+1
    }
    #Add the country with smallest SSE to the group
    choice=Pool[[Index]]
    Group=Groupdemogdata_V2(Pool[[which(names(Pool)==target$label)]],
                         choice,choice$year[1],choice$year[length(choice$year)])
    Best[[component_size,2]]=Pool[which(names(Pool) %in% Group$label)]
    LOSS[component_size,2]=error
    #Update the candidate pool
    PoolStatus[Index]=0
    size=size+1
    ## DSA till size reaches CoG
    while(size<COG)
    {
      message("Current groupsize is ",size)
      Signal=0 #if Adjustment Signal=1, omit subsequent codes and go back to this line
      ## Deletion
      ## when size>3, judge whether we could delete one country to obtain better results for (size-1) group
      if(size>2 && size<COG)
      {
        error=Inf
        Current_Group=Best[[component_size,size]]
        for(country in Current_Group)
        {
          label=country$label
          if(label!=target$label)
          {
            temp=Current_Group[names(Current_Group)!=label]
            group=Pool[[which(names(Pool)==target$label)]]#target
            for( member in temp)
            {
              if(member$label!=target$label)
              {
                group=Groupdemogdata_V2(group,member,member$year[1],member$year[length(member$year)])
              }
            }
            error_new=loss(prefit_v2(group,Target=target,component_size,Pool))[[LossFunction]]
            if(error_new<error)
            {
              error=error_new
              Choice=label
            }
          }
        }
        if(error<LOSS[component_size,size-1])
          ## when deletion can yield better SSE for size-1 group
        {
          Group=Current_Group[which(names(Current_Group)!=Choice)]
          Best[[component_size,size-1]]=Pool[which(names(Pool) %in% names(Group))]
          ## Waring message
          if(length(Best[[component_size,size-1]])!=(size-1))
          {
            message(Choice)
            message("current size= ",size," but deleted group size= ",length(Best[[component_size,size-1]]))
          }
          LOSS[component_size,size-1]=error
          ## the deleted country would go back to the candidate pool
          index=which(names(Pool)==Choice)
          PoolStatus[index]=1
          size=size-1
          Signal=1
        }
      }

      #########################################################################
      ## Substitution
      ## When deletion is not possible(Signal=0), judge whether we could substitute one country to obtain better results for (size) group
      if(Signal==0)
      {
        error=Inf
        Current_Group=Best[[component_size,size]]
        ## first delete one country from original group
        for(country1 in Current_Group)
        {
          label1=country1$label
          if(label1!=target$label)
          {
            temp=Current_Group[names(Current_Group)!=label1]
            group=Pool[[which(names(Pool)==target$label)]]#target
            for( member in temp)
            {
              if(member$label!=target$label)
              {
                group=Groupdemogdata_V2(group,member,member$year[1],member$year[length(member$year)])
              }
            }
            for (country2 in Pool)
            {
              label2=country2$label
              index2=which(names(Pool)==label2)
              if(PoolStatus[index2]==1)
              {
                group_new=Groupdemogdata_V2(group,country2,country2$year[1],country2$year[length(country2$year)])
                error_new=loss(prefit_v2(group_new,Target=target,component_size,Pool))[[LossFunction]]
                if(error_new<error)
                {
                  error=error_new
                  Choice1=label1
                  Choice2=label2
                }
              }
            }
          }
        }
        if(error<LOSS[component_size,size])
          ## when substitution can achieve smaller sse for current size group
        {
          ## in Best[component_size,size], substitute Choice1 with Choice2
          ## in Pool, delete Choice2 and add Choice1
          ## Update Best[[component_size,size]]
          message("Current size is ",size," Choice1 is ",Choice1," Choice2 is ",Choice2)
          Index=which(names(Current_Group)==Choice1)
          Best[[component_size,size]][[Index]]=Pool[[which(names(Pool)==Choice2)]]
          names(Best[[component_size,size]])=sapply(Best[[component_size,size]],GetLabel)
          message(names(Best[[component_size,size]]))
          ## Update Pool
          Index=which(names(Pool)==Choice1)
          PoolStatus[Index]=1
          Index=which(names(Pool)==Choice2)
          PoolStatus[Index]=0
          LOSS[component_size,size]=error
          Signal=1
        }
      }
      #########################################################################
      ## Addition (if Signal=0)
      ## based on current best (size) group, judge whether we could add one country to obtain (size+1) group
      if(Signal==0)
      {
        ## Obrain the grouped demogdata of Best[component_size,size]
        Group=Pool[[which(names(Pool)==target$label)]]#target
        for(country in Best[[component_size,size]])
        {
          label=country$label
          if(label!=target$label)
            Group=Groupdemogdata_V2(Group,country,country$year[1],country$year[length(country$year)])
        }
        ## Define addition function, calculate SSE if we add-in one country
        addition<-function(Obj)
        {
          index=which(names(Pool)==Obj$label)
          if(PoolStatus[index]==0)
          {
            Out=Inf
          }
          else
          {
            temp=Groupdemogdata_V2(Group,Obj,Obj$year[1],Obj$year[length(Obj$year)])
            Out=loss(prefit_v2(temp,Target=target,component_size,Pool))[[LossFunction]]
          }
          return(Out)
        }
        AddResult=sapply(Pool,addition)
        Index=which(AddResult==min(AddResult))
        #Add the country with smallest SSE to the group
        choice=Pool[[Index[1]]]
        Group_label=Groupdemogdata_V2(Group,choice,choice$year[1],choice$year[length(choice$year)])$label
        Best[[component_size,size+1]]=Pool[which(names(Pool) %in% Group_label)]
        LOSS[component_size,size+1]=min(AddResult)
        #Update the candidate pool
        PoolStatus[Index]=0
        size=size+1
      }
    }
  }
  Out=list(target=target,COG=COG,N=N,LossFunction=LossFunction,Poolnames=Poolnames,Best=Best,LOSS=LOSS)
  return(Out)
}
