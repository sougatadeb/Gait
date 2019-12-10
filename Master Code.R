clVar <- function() {
  env <- parent.frame()
  rm(list = setdiff( ls(all.names=TRUE, env = env), lsf.str(all.names=TRUE, env = env)),envir = env)
}
clVar()
require(zoo)
warp <- function(str2,newRow){
  str3 <- rep(0,newRow*ncol(str2))
  dim(str3) <- c(newRow,ncol(str2))
  incr <- (nrow(str2)-1)/(newRow-1)
  for (strInd in 1:newRow) {
    useInd <- 1 + (strInd-1)*incr
    str3[strInd,1] <- str2[floor(useInd),1]*(1-useInd+floor(useInd)) +
      str2[min(nrow(str2),floor(useInd+incr)),1]*(useInd-floor(useInd))
    str3[strInd,2] <- str2[floor(useInd),2]*(1-useInd+floor(useInd)) +
      str2[min(nrow(str2),floor(useInd+incr)),2]*(useInd-floor(useInd))
    str3[strInd,3] <- str2[floor(useInd),3]*(1-useInd+floor(useInd)) +
      str2[min(nrow(str2),floor(useInd+incr)),3]*(useInd-floor(useInd))
  }
  return(str3)
}
RoR.Calc <- function(prStr,wts) {
  cc <- cbind((prStr[,4]-mean(prStr[,4]))^2*wts[1]+(prStr[,5]-mean(prStr[,5]))^2*wts[2]+(prStr[,6]-mean(prStr[,6]))^2*wts[3],
              (prStr[,1]-prStr[,4])^2*wts[1]+(prStr[,2]-prStr[,5])^2*wts[2]+(prStr[,3]-prStr[,6])^2*wts[3],rep(1,nrow(prStr)))
  dd <- cc[order(cc[,2]),]
  ee <- apply(dd,2,cumsum)
  return (max(0.1,max((1-ee[,2,drop=FALSE]/(ee[,1,drop=FALSE]+0.000001))*sqrt(ee[,3,drop=FALSE]/nrow(prStr)))))
}
run.Train <- function(Typ, rollP, RRCut, nPoint, nFreq){
  trainAll <- read.csv(paste("Clean/",Typ,"_TRN.csv",sep = ""), header=T, na.strings=c(""," ","NA"))
  userList <- unique(trainAll[,4])
  prdUsr <- rep(0,length(userList))
  trainCln <- NULL
  for (user in 1:length(userList)){
    trainR <- trainAll[which(trainAll[,4]==userList[user]),]
    trainR <- trainR[-c(1:(3*nFreq)),]
    if (rollP > 1) {
      train <- data.frame(cbind(rollapply(trainR[,1],rollP,mean),
                                rollapply(trainR[,2],rollP,mean),
                                rollapply(trainR[,3],rollP,mean)))
    } else {
      train <- trainR
    }
    rm(trainR)
    colnames(train) <- c("X","Y","Z")
    trainFin <- NULL
    acfX <- acf(train[1:(3*nFreq),1],lag.max = round(1.8*nFreq), plot = F)
    acfY <- acf(train[1:(3*nFreq),2],lag.max = round(1.8*nFreq), plot = F)
    acfZ <- acf(train[1:(3*nFreq),3],lag.max = round(1.8*nFreq), plot = F)
    acfA <- apply(cbind(acfX$acf,acfY$acf,acfZ$acf),1,mean)
    prdUsr[user] <- floor(nFreq*1/2)+which.max(acfA[(floor(nFreq*1/2)+1):(round(1.8*nFreq))])
    testIt <- train[   1:prdUsr[user] ,1:3]
    train  <- train[-c(1:prdUsr[user]),1:3]
    trainFin <- cbind(warp(testIt,nPoint),rep(nrow(testIt),nPoint))
    colnames(trainFin) <- c("X1","X2","X3","X7")
    while (nrow(train) >= 3*nFreq) {
      acfX <- acf(train[1:(3*nFreq),1],lag.max = round(1.8*nFreq), plot = F)
      acfY <- acf(train[1:(3*nFreq),2],lag.max = round(1.8*nFreq), plot = F)
      acfZ <- acf(train[1:(3*nFreq),3],lag.max = round(1.8*nFreq), plot = F)
      acfA <- apply(cbind(acfX$acf,acfY$acf,acfZ$acf),1,mean)
      prdUsr[user] <- floor(nFreq*1/2)+which.max(acfA[(floor(nFreq*1/2)+1):(round(1.8*nFreq))])
      testsWa <- warp(testIt,prdUsr[user])
      trainWa <- rbind(train[1:prdUsr[user],1:3],train[1:prdUsr[user],1:3])
      ccfX <- ccf(trainWa[,1],testsWa[,1], lag.max = prdUsr[user], plot = F)
      ccfY <- ccf(trainWa[,2],testsWa[,2], lag.max = prdUsr[user], plot = F)
      ccfZ <- ccf(trainWa[,3],testsWa[,3], lag.max = prdUsr[user], plot = F)
      ccfA <- apply(cbind(ccfX$acf,ccfY$acf,ccfZ$acf),1,mean)
      if (max(ccfA) > RRCut) {
        continue <- ifelse(which.max(ccfA)>=prdUsr[user],which.max(ccfA)-prdUsr[user]+1,which.max(ccfA)+1)
        testIt <- trainWa[(continue+0):(continue+prdUsr[user]-1),1:3,drop=F]
        trainFin <- rbind(cbind(warp(testIt,nPoint),rep(prdUsr[user],nPoint)),trainFin)
      } 
      train <- train[-c(1:prdUsr[user]),1:3]
    }
    trainUsr <- cbind(trainFin,rep(userList[user],nrow(trainFin)))
    rm(trainFin)
    colnames(trainUsr) <- c("X","Y","Z","Period","User")
    trainCln <- rbind(trainCln,trainUsr)
    rm(trainUsr)
  }
  write.csv(trainCln,paste("Clean/",Typ,"_CLNAVG.csv",sep = ""), row.names = F)
}
run.Test <- function(Typ, Repli, rollP, RRCut, nPoint, nFreq) {
  trainCln <- read.csv(paste("Clean/",Typ,"_CLNAVG.csv",sep = ""), header=T, na.strings=c(""," ","NA"))
  testsAll <- read.csv(paste("Clean/",Typ,"_TST.csv",sep = ""), header=T, na.strings=c(""," ","NA"))
  steps <- floor(nrow(trainCln)/nPoint)
  userList <- unique(trainCln[,5])
  for (iter in 1:Repli) {
    dummy <- 0
    while (dummy == 0) {
      dPoint <- sample(1:(nrow(testsAll)-4*nFreq-rollP+1),1)
      if (testsAll[dPoint,4] == testsAll[(dPoint+4*nFreq+rollP-1),4]) {
        dummy <- 1
      }      
    }
    testsR <- testsAll[dPoint:(dPoint+4*nFreq+rollP-1),1:3,drop=F]
    if (rollP > 1) {
      testRo <- data.frame(cbind(rollapply(testsR[,1],rollP,mean),
                                 rollapply(testsR[,2],rollP,mean),
                                 rollapply(testsR[,3],rollP,mean)))
    } else {
      testRo <- testsR
    }
    rm(testsR)
    colnames(testRo) <- c("X","Y","Z")
    acfX <- acf(testRo[,1],lag.max = round(1.8*nFreq), plot = F)
    acfY <- acf(testRo[,2],lag.max = round(1.8*nFreq), plot = F)
    acfZ <- acf(testRo[,3],lag.max = round(1.8*nFreq), plot = F)
    acfA <- apply(cbind(acfX$acf,acfY$acf,acfZ$acf),1,mean)
    prdTst <- floor(nFreq*1/2)+which.max(acfA[(floor(nFreq*1/2)+1):(round(1.8*nFreq))])
    testIt <- warp(testRo[1:prdTst,1:3,drop=F],nPoint)
    testI2 <- warp(testRo[(prdTst+1):(2*prdTst),1:3,drop=F],nPoint)
    
    RR <- rep(0,steps)
    UR <- rep(0,steps)
    userC <- 0
    userV <- rep(0,length(userList))
    userM <- rep(0,length(userList))
    user3 <- rep(0,length(userList))
    for(j in 1:steps){
      trainWa <- rbind(trainCln[((j-1)*nPoint+1):(j*nPoint),1:3, drop=F],trainCln[((j-1)*nPoint+1):(j*nPoint),1:3, drop=F])
      ccfX <- ccf(trainWa[,1],testIt[,1], lag.max = nPoint, plot = F)
      ccfY <- ccf(trainWa[,2],testIt[,2], lag.max = nPoint, plot = F)
      ccfZ <- ccf(trainWa[,3],testIt[,3], lag.max = nPoint, plot = F)
      ccfA <- apply(cbind(ccfX$acf,ccfY$acf,ccfZ$acf),1,mean)
      ccfX <- ccf(trainWa[,1],testI2[,1], lag.max = nPoint, plot = F)
      ccfY <- ccf(trainWa[,2],testI2[,2], lag.max = nPoint, plot = F)
      ccfZ <- ccf(trainWa[,3],testI2[,3], lag.max = nPoint, plot = F)
      ccfB <- apply(cbind(ccfX$acf,ccfY$acf,ccfZ$acf),1,mean)
      if (max(ccfA) > RRCut) {
        continue <- ifelse(which.max(ccfA)>=nPoint,which.max(ccfA)-nPoint+1,which.max(ccfA)+1)
        paired <- cbind(trainWa[(continue+0):(continue+nPoint-1),],testIt)
        X <- as.matrix(cbind(rep(1,nrow(paired)),paired[,1:3]))
        Y <- as.matrix(paired[,4:6])
        b <- solve(crossprod(X))%*%crossprod(X,Y)
        ROR1 <- sqrt(RoR.Calc(cbind(X%*%b,paired[,4:6]),c(1,1,1))*RoR.Calc(paired,c(1,1,1)))
        continue <- ifelse(which.max(ccfB)>=nPoint,which.max(ccfB)-nPoint+1,which.max(ccfB)+1)
        paired <- cbind(trainWa[(continue+0):(continue+nPoint-1),],testI2)
        X <- as.matrix(cbind(rep(1,nrow(paired)),paired[,1:3]))
        Y <- as.matrix(paired[,4:6])
        b <- solve(crossprod(X))%*%crossprod(X,Y)
        ROR2 <- sqrt(RoR.Calc(cbind(X%*%b,paired[,4:6]),c(1,1,1))*RoR.Calc(paired,c(1,1,1)))
        RR[j] <- max(ifelse(is.na(ROR1)==T,0,ROR1),ifelse(is.na(ROR2)==T,0,ROR2))
      } else {
        RR[j] <- 0
      }
      UR[j] <- trainCln[((j-1)*nPoint+1),5]
      if (j==1) {
        UserP <- RR[j]
      } else if (UR[j]==UR[(j-1)]) {
        UserP <- c(UserP,RR[j])
        if (j == steps) {
          userC <- userC + 1
          userV[userC] <- mean(UserP[which(UserP>=quantile(UserP,probs = 0.9))])
          user3[userC] <- mean(UserP[order(UserP, decreasing=TRUE)[1:3]])
          userM[userC] <- max(UserP)
        }
      } else {
        userC <- userC + 1
        userV[userC] <- mean(UserP[which(UserP>=quantile(UserP,probs = 0.9))])
        user3[userC] <- mean(UserP[order(UserP, decreasing=TRUE)[1:3]])
        userM[userC] <- max(UserP)
        UserP <- RR[j]
      }
    }
    bestRSQ <- userM[order(userM, decreasing=TRUE)[1:2]]
    bestUSR <- userList[order(userM, decreasing=TRUE)[1:2]]
    bestP90 <- userList[which.max(userV)]
    bestR03 <- userList[which.max(user3)]
    URT <- UR[which(RR>=quantile(RR,probs = 0.95))]
    RRT <- RR[which(RR>=quantile(RR,probs = 0.95))]
    aa <- as.matrix(data.frame(table(URT)))
    bestTOP <- aa[which.max(aa[,2]),1]
    UserID <- c(testsAll[dPoint,4],bestP90,bestTOP,bestUSR,bestR03,round(max(userV),4),
                round(max(RRT[which(URT==bestTOP)]),4),round(bestRSQ,4),round(max(user3),4),prdTst,dPoint,Typ)
    dim(UserID) <- c(1,14)
    colnames(UserID) <- c("Act","bestP90","bestTOP","bestMax","best2ND","bestOf3",
                          "RSQP90","RSQTOP","RSQMax","RSQ2ND","RSQ3","Period","Sample","Activity")
    UserID <- cbind(UserID,round(t(userM),3))
    if (iter == 1) {
      finDat <- UserID
    } else {
      finDat <- rbind(finDat, UserID)
    }
    if (iter%%5 == 0) {
      write.csv(finDat,paste("Results/",Typ,"_",Repli,"_",rollP,"_","NEWAVG.csv",sep = ""), row.names = F)
    }
  }
  return(finDat)
}

    
setwd("C:/1 Course/KE5006. AppRes (Done)/1. Data/1. MAREA/")
run.Train(Typ = 'IFSRL', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'IFSRR', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'IFSRT', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'IFSRW', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'IFSWL', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'IFSWR', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'IFSWT', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'IFSWW', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'OSSRL', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'OSSRR', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'OSSRT', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'OSSRW', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'OSSWL', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'OSSWR', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'OSSWT', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'OSSWW', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'TRFRL', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'TRFRR', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'TRFRT', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'TRFRW', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'TRFWL', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'TRFWR', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'TRFWT', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'TRFWW', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'TRSWL', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'TRSWR', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'TRSWT', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)
run.Train(Typ = 'TRSWW', rollP = 5, RRCut = 0.00, nPoint = 100, nFreq = 128)

setwd("C:/1 Course/KE5006. AppRes (Done)/1. Data/2. HAPT/")
run.Train(Typ = 'WALKDN', rollP = 3, RRCut = 0.00, nPoint = 40, nFreq = 50)
run.Train(Typ = 'WALKIT', rollP = 3, RRCut = 0.00, nPoint = 40, nFreq = 50)
run.Train(Typ = 'WALKUP', rollP = 3, RRCut = 0.00, nPoint = 40, nFreq = 50)

setwd("C:/1 Course/KE5006. AppRes (Done)/1. Data/1. MAREA/")
IFSRL <- run.Test(Typ = 'IFSRL', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
IFSRR <- run.Test(Typ = 'IFSRR', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
IFSRT <- run.Test(Typ = 'IFSRT', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
IFSRW <- run.Test(Typ = 'IFSRW', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
IFSWL <- run.Test(Typ = 'IFSWL', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
IFSWR <- run.Test(Typ = 'IFSWR', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
IFSWT <- run.Test(Typ = 'IFSWT', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
IFSWW <- run.Test(Typ = 'IFSWW', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
OSSRL <- run.Test(Typ = 'OSSRL', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
OSSRR <- run.Test(Typ = 'OSSRR', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
OSSRT <- run.Test(Typ = 'OSSRT', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
OSSRW <- run.Test(Typ = 'OSSRW', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
OSSWL <- run.Test(Typ = 'OSSWL', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
OSSWR <- run.Test(Typ = 'OSSWR', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
OSSWT <- run.Test(Typ = 'OSSWT', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
OSSWW <- run.Test(Typ = 'OSSWW', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
TRFRL <- run.Test(Typ = 'TRFRL', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
TRFRR <- run.Test(Typ = 'TRFRR', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
TRFRT <- run.Test(Typ = 'TRFRT', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
TRFRW <- run.Test(Typ = 'TRFRW', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
TRFWL <- run.Test(Typ = 'TRFWL', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
TRFWR <- run.Test(Typ = 'TRFWR', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
TRFWT <- run.Test(Typ = 'TRFWT', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
TRFWW <- run.Test(Typ = 'TRFWW', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
TRSWL <- run.Test(Typ = 'TRSWL', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
TRSWR <- run.Test(Typ = 'TRSWR', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
TRSWT <- run.Test(Typ = 'TRSWT', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)
TRSWW <- run.Test(Typ = 'TRSWW', Repli = 1000, rollP = 5, RRCut = 0.0, nPoint = 100, nFreq = 128)

setwd("C:/1 Course/KE5006. AppRes (Done)/1. Data/2. HAPT/")
WALKDN <- run.Test(Typ = 'WALKDN', Repli = 3000, rollP = 3, RRCut = 0.0, nPoint = 40, nFreq = 50)
WALKIT <- run.Test(Typ = 'WALKIT', Repli = 3000, rollP = 3, RRCut = 0.0, nPoint = 40, nFreq = 50)
WALKUP <- run.Test(Typ = 'WALKUP', Repli = 3000, rollP = 3, RRCut = 0.0, nPoint = 40, nFreq = 50)

setwd("C:/1 Course/KE5006. AppRes/1. Data (Done)/1. MAREA/")
COMBINED1 <- rbind(IFSRL,IFSRR,IFSRW,IFSWL,IFSWR,IFSWW,TRFRL,TRFRR,TRFRW,TRFWL,TRFWR,TRFWW,TRSWL,TRSWR,TRSWW)
write.csv(COMBINED1,"Results/COMBINED1_NEW.csv", row.names = F)
COMBINED4 <- rbind(OSSRL,OSSRR,OSSRW,OSSWL,OSSWR,OSSWW,OSSRT,OSSWT)
write.csv(COMBINED4,"Results/COMBINED4_NEW.csv", row.names = F)
COMBINED2 <- rbind(IFSRT,IFSWT,TRFRT,TRFWT,TRSWT)
write.csv(COMBINED2,"Results/COMBINED2_NEW.csv", row.names = F)
COMBINED3 <- rbind(WALKDN,WALKIT,WALKUP)
write.csv(COMBINED3,"Results/COMBINED3_NEW.csv", row.names = F)

