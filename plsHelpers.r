mapPCAPLS<-function(newData, roplsObject)
{
  # This fucntion calculate scores (T) transformation of new data. So basically the new dataset will be mapped to the old model space
  
  # check the all the arguments is present
  if(missing(newData) || missing(roplsObject))
  {
    stop("All the argument (newData: new data matrix, roplsObject: ropls object) should be present!")
  }
  if(roplsObject@typeC=="OPLS-DA")
  {
    cat("We are assuming oplsda input now!")
    # Recalculate orthogonal scores
    Tortho<-roplsObject@suppLs$xModelMN %*% roplsObject@orthoWeightMN %*% solve(t(roplsObject@orthoWeightMN) %*% roplsObject@orthoWeightMN)
    # Recalculate orthogonal loading
    Portho<- t(roplsObject@suppLs$xModelMN) %*% Tortho %*% solve(t(Tortho) %*% Tortho) 
    # reconstruct the data based on the X=TP'
    # and substract the data fronm the original matrix
    xtrainCheck <- roplsObject@suppLs$xModelMN - Tortho %*% t(Portho)
    # check if we get the same results compared to the original roplsda
    # we might encounter some very minor differences (some decimal points)
    # so we use round!
    calculationCheck<-all(round(roplsObject@scoreMN,5)==round(xtrainCheck%*%roplsObject@weightMN,5))
    if(!calculationCheck)
    {
      stop("We failed to calculate back the predictive scores of the original oplsda! There is something wrong!")
    }
    # Scale the new data based on the center and sd of the old data
    X<-scale(newData,center = roplsObject@xMeanVn,roplsObject@xSdVn)
    # calculate orthogonal scores
    Tortho<-X %*% roplsObject@orthoWeightMN %*% solve(t(roplsObject@orthoWeightMN) %*% roplsObject@orthoWeightMN)
    # reconstruct the data based on the X=TP'
    # and substract the data fronm the original matrix
    xtest <- X - Tortho %*% t(Portho)
    # calculate the predictive scores (T) for the new data
    xtestScoreMN<-xtest %*% roplsObject@weightMN
    # return the results
    return(xtestScoreMN)
  }else if (roplsObject@typeC=="PCA")
  {
    cat("We are assuming PCA input now!")
    
    # reconstruct the original PCA
    calculationCheck<-all(round(roplsObject@scoreMN,5)==round(roplsObject@suppLs$xModelMN %*% roplsObject@loadingMN,5))
    
    if(!calculationCheck)
    {
      stop("We failed to calculate back the  scores of the original PCA! There is something wrong!")
    }
    # Scale the new data based on the center and sd of the old data
    X<-scale(newData,center = roplsObject@xMeanVn,roplsObject@xSdVn)
    # calculate the scores (T) for the new data
    xtestScoreMN<-X %*% roplsObject@loadingMN
    # return the results
    return(xtestScoreMN)
    
  } else if(roplsObject@typeC=="PLS-DA")
  {
    
    cat("We are assuming PLS-DA input now!")
    
    # reconstruct the original PLSDA
    calculationCheck<-all(round(roplsObject@scoreMN,5)==round(roplsObject@suppLs$xModelMN %*% roplsObject@weightStarMN,5))
    
    if(!calculationCheck)
    {
      stop("We failed to calculate back the  scores of the original PCA! There is something wrong!")
    }
    # Scale the new data based on the center and sd of the old data
    X<-scale(newData,center = roplsObject@xMeanVn,roplsObject@xSdVn)
    # calculate the scores (T) for the new data
    xtestScoreMN<-X %*% roplsObject@weightStarMN
    # return the results
    return(xtestScoreMN)
    
  }

}

