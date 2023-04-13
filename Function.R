# Packages ----
my_packages <- c("tidyverse", "data.table",
                 "GSEABase",
                 "survival", "survminer",
                 "car",
                 "My.stepwise",
                 "randomForest","glmnet",
                 "plyr", "ROCR",
                 "rsample", "withr",
                 "R.utils",
                 "snowfall", "parallel")
pacman::p_load(char = my_packages)

# basic----
basicFilter <- function(test_exp){
  passed_genes <- apply(test_exp, 1, function(genes){
    data.frame(
      mean = mean(genes),
      sd = sd(genes)
    )
  })
  passed_genes <- do.call(rbind, passed_genes)
  passed_genes <- passed_genes[passed_genes$mean > 1 &
                                 passed_genes$sd > 0.1,]
  test_exp2 <- test_exp[rownames(test_exp) %in% rownames(passed_genes),] 
  return(test_exp2)
}

tumorFilter <- function(test_exp, test_survival){
  group = as.numeric(substr(colnames(test_exp), nchar(colnames(test_exp))-2, nchar(colnames(test_exp))-1)) < 10
  test_exp2 <- test_exp[,group]
  test_survival2 <- test_survival[match(colnames(test_exp2), test_survival$sample),]
  return(list(test_exp2,
              test_survival2))
}

# wilcox----
wilcoxDEGs_fast <- function(exp_df, clinical1){
  # exp_df in log2 format
  exp_df <- 2**(exp_df) - 1
  sfInit(parallel = TRUE, cpus = 6)
  sfExport("exp_df", "clinical1")  
  re <- sfLapply(1:nrow(exp_df), function(x){
    id1 <- clinical1[clinical1$group %in% levels(clinical1$group)[1],]$patient
    id2 <- clinical1[clinical1$group %in% levels(clinical1$group)[2],]$patient
    group1 <- exp_df[,colnames(exp_df) %in% id1]
    group2 <- exp_df[,colnames(exp_df) %in% id2]
    p <- wilcox.test(as.numeric(group2[x,]), as.numeric(group1[x,]))$p.value
    logFC <- log2(mean(as.numeric(group2[x,])) / mean(as.numeric(group1[x,])))
    c(p, logFC)
  })
  re <- as.data.frame(do.call(rbind, re))
  rownames(re) <- rownames(exp_df)
  colnames(re) <- c("P.value", "LogFC")
  re$FDR <- p.adjust(re$P.value, method = "fdr") 
  sfStop()
  return(re)
}

wilcoxDEGs <- function(exp_df, clinical1){
  # exp_df in log2 format
  exp_df <- 2**(exp_df) - 1
  re <- lapply(1:nrow(exp_df), function(x){
    id1 <- clinical1[clinical1$group %in% levels(clinical1$group)[1],]$patient
    id2 <- clinical1[clinical1$group %in% levels(clinical1$group)[2],]$patient
    group1 <- exp_df[,colnames(exp_df) %in% id1]
    group2 <- exp_df[,colnames(exp_df) %in% id2]
    p <- wilcox.test(as.numeric(group2[x,]), as.numeric(group1[x,]))$p.value
    logFC <- log2(mean(as.numeric(group2[x,])) / mean(as.numeric(group1[x,])))
    c(p, logFC)
  })
  re <- as.data.frame(do.call(rbind, re))
  rownames(re) <- rownames(exp_df)
  colnames(re) <- c("P.value", "LogFC")
  re$FDR <- p.adjust(re$P.value, method = "fdr") 
  return(re)
}

# Cox----
coxFilter_fast <- function(test_exp3, test_survival3){
  sfInit(parallel = TRUE, cpus = 6)
  test_surv <- with(test_survival3, Surv(OS.time, OS))
  sfExport("test_exp3", "test_survival3") 
  sfLibrary(survival)
  re <- sfApply(test_exp3, 1, function(gene){
    res.cox <- summary(coxph(test_surv ~ gene))
    # res.cox$waldtest[3]
    return(res.cox$waldtest[3])
  })
  sfStop()
  return(re)
}

coxFilter <- function(test_exp3, test_survival3){
  test_surv <- with(test_survival3, Surv(OS.time, OS))
  re <- apply(test_exp3, 1, function(gene){
    res.cox <- summary(coxph(test_surv ~ gene))
    # res.cox$waldtest[3]
    return(res.cox$waldtest[3])
  })
  return(re)
}

# Model Construction----
calAUC <- function(re){
  pred <- prediction(re[,1], re[,2]) # predict events by prob_min
  auc <- performance(pred,"auc")@y.values[[1]] # get auc value
  return(auc)
}

myLassoLike <- function(test_exp5, 
                        test_survival4,
                        seeds = c(1,2),
                        alphas = c(0,0.5,1)){
  re <- lapply(seeds, function(seed){
    with_seed(seed = seed, svfold_cv <- vfold_cv(test_survival4, v = 3, repeats = 5,                                               strata = OS))
    svfold_cv <- svfold_cv[[1]]
    a <- 1
    lapply(svfold_cv, function(cv){
      # cv = svfold_cv[[2]]
      # Split data 
      train_survival <- training(cv)
      val_survival <- testing(cv)
      train_exp <- test_exp5[,match(train_survival$sample, colnames(test_exp5))]
      val_exp <- test_exp5[,match(val_survival$sample, colnames(test_exp5))]
      # model construction
      auc1 <- lapply(alphas, function(alpha){
        with_seed(seed = 66, lasso_CV <- cv.glmnet(x=t(train_exp), y=train_survival$OS, nlambda = 1000,alpha = alpha))
        ## Train
        lasso.pre.train <- as.data.frame(predict(lasso_CV, newx=t(train_exp), s=lasso_CV$lambda.min))
        lasso.pre.train$True <- train_survival$OS
        train.auc <- tryCatch(calAUC(lasso.pre.train),
                              error = function(x){0})
        ## Test
        lasso.pre.test <- as.data.frame(predict(lasso_CV, newx=t(val_exp), s=lasso_CV$lambda.min))
        lasso.pre.test$True <- val_survival$OS
        test.auc <- tryCatch(calAUC(lasso.pre.test),
                             error = function(x){0})
        index <- a
        # collect genes
        best_model <- glmnet(x=t(train_exp), y=train_survival$OS, alpha = alpha, lambda = lasso_CV$lambda.min)
        best_model_re <- as.matrix(coef(best_model))
        best_model_re <- best_model_re[,1]
        best_model_re <- best_model_re[-1]
        selected_lasso_genes <- names(best_model_re[best_model_re!=0])
        selected_lasso_genes <- paste0(selected_lasso_genes, collapse = ", ")
        # print(alpha)
        c(train.auc,
          test.auc,
          seed, index, alpha, selected_lasso_genes)
      })
      # print(a)
      a <<- a + 1
      auc1 <- do.call(rbind, auc1)
      return(auc1)
    })
  })
  return(re)
}


myCoxLike <- function(test_exp5, 
                      test_survival4,
                      seeds = c(1,2)){
  re <- lapply(seeds, function(seed){
    # test_survival4 = test_survival5
    with_seed(seed = seed, svfold_cv <- vfold_cv(test_survival4, v = 3, repeats = 5, strata = OS))
    svfold_cv <- svfold_cv[[1]]
    a <- 1
    re <- lapply(svfold_cv, function(cv){
      # cv = svfold_cv[[1]]
      # Split data 
      train_survival <- training(cv)
      val_survival <- testing(cv)
      train_exp <- test_exp5[,match(train_survival$sample, colnames(test_exp5))]
      val_exp <- test_exp5[,match(val_survival$sample, colnames(test_exp5))]
      # model construction
      variables <- rownames(train_exp)
      my_formula <- formula(paste0("Surv(OS.time, OS) ~ ", paste0(variables, collapse = " + ")))
      coxph_data <- cbind(train_survival, t(train_exp))
      cox_model <- coxph(formula = my_formula, data = coxph_data, method = "efron")
      # filter cox for accelerating step.wise
      cox_re <- coxFilter(train_exp, train_survival)
      cox_re2 <- cox_re[cox_re < .1]
      p_val <- .05
      while(length(cox_re2) > 15) {
        cox_re2 <- cox_re2[cox_re2 < p_val]
        p_val = p_val / 2
      }
      coxph_data2 <- coxph_data[,colnames(coxph_data) %in% c(names(cox_re2), colnames(coxph_data)[1:5])]
      # coxph_data2 <- cbind(coxph_data[,1:5], coxph_data2)
      stepwise_cox_model <- tryCatch({withTimeout(stepwise_cox_model <- stepwiseCox(Time = "OS.time", Status = "OS",
                                                                                    variable.list = colnames(coxph_data2)[-c(1:5)], data = coxph_data2), timeout = 80)}, error = function(x){return(cox_model)})
      ## Train
      cox.pre.train <- as.data.frame(predict(cox_model, data = t(train_exp)))
      cox.pre.train$True <- train_survival$OS
      train.auc <- tryCatch(calAUC(cox.pre.train),
                            error = function(x){0})
      ## Test
      val_meta <- cbind(as.data.frame(t(val_exp)),
                        val_survival[,4:5]) 
      cox.pre.test <- as.data.frame(predict(cox_model, newdata = val_meta))
      cox.pre.test$True <- val_survival$OS
      test.auc <- tryCatch(calAUC(cox.pre.test),
                           error = function(x){0})
      index <- a
      # print(index)
      a <<- a + 1
      sel_genes <- names(cox_model$coefficients)
      sel_genes <- paste0(sel_genes, collapse = ", ")
      re1 <- c(train.auc,
               test.auc,
               seed, index, "raw_cox", sel_genes)
      ## Train
      cox.pre.train <- as.data.frame(predict(stepwise_cox_model, data = t(train_exp)))
      cox.pre.train$True <- train_survival$OS
      train.auc <- tryCatch(calAUC(cox.pre.train),
                            error = function(x){0})
      ## Test
      cox.pre.test <- as.data.frame(predict(stepwise_cox_model, newdata =  val_meta))
      cox.pre.test$True <- val_survival$OS
      test.auc <- tryCatch(calAUC(cox.pre.test),
                           error = function(x){0})
      sel_genes2 <- names(stepwise_cox_model$coefficients)
      sel_genes2 <- paste0(sel_genes2, collapse = ", ")                           
      re2 <- c(train.auc,
               test.auc,
               seed, index, "stepwise_cox", sel_genes2)
      data.frame(
        rbind(re1, re2)
      )
    })
  })
  return(re)
}

stepwiseCox <- function (Time = NULL, T1 = NULL, T2 = NULL, Status = NULL, variable.list, 
                         in.variable = "NULL", data, sle = 0.15, sls = 0.15, vif.threshold = 999) 
{
  univar.pvalue <- NULL
  temp.model <- NULL
  if (is.null(T2)) {
    initial.model <- coxph(as.formula(paste("Surv(", Time, 
                                            ", ", Status, ") ~ ", paste(in.variable, collapse = "+"), 
                                            sep = "")), data = data, method = "efron")
  }
  else if (is.null(T2)) {
    initial.model <- coxph(as.formula(paste("Surv(", Time, 
                                            ", ", Status, ") ~ ", paste(in.variable, collapse = "+"), 
                                            sep = "")), data = data, method = "efron")
  }
  else if (is.null(Time)) {
    initial.model <- coxph(as.formula(paste("Surv(", T1, 
                                            ", ", T2, ", ", Status, ") ~ ", paste(in.variable, 
                                                                                  collapse = "+"), sep = "")), data = data, method = "efron")
  }
  else if (is.null(Time)) {
    initial.model <- coxph(as.formula(paste("Surv(", T1, 
                                            ", ", T2, ", ", Status, ") ~ ", paste(in.variable, 
                                                                                  collapse = "+"), sep = "")), data = data, method = "efron")
  }
  if (is.null(initial.model$coefficients)) {
    for (i in 1:length(variable.list)) {
      if (is.null(T2)) {
        uni.model <- coxph(as.formula(paste("Surv(", 
                                            Time, ", ", Status, ") ~ ", variable.list[i], 
                                            sep = "")), data = data, method = "efron")
      }
      if (is.null(Time)) {
        uni.model <- coxph(as.formula(paste("Surv(", 
                                            T1, ", ", T2, ", ", Status, ") ~ ", variable.list[i], 
                                            sep = "")), data = data, method = "efron")
      }
      univar.pvalue[i] <- summary(uni.model)$coefficients[5]
    }
    variable.list1 <- variable.list[univar.pvalue <= 0.9 & 
                                      !is.na(univar.pvalue)]
    univar.pvalue1 <- univar.pvalue[univar.pvalue <= 0.9 & 
                                      !is.na(univar.pvalue)]
    uni.x <- variable.list1[which.min(univar.pvalue1)]
    if (length(uni.x) > 0) {
      if (is.null(T2)) {
        formula <- as.formula(paste("Surv(", Time, ", ", 
                                    Status, ") ~ ", uni.x, sep = ""))
        temp.model <- coxph(formula, data = data, method = "efron")
      }
      if (is.null(Time)) {
        formula <- as.formula(paste("Surv(", T1, ", ", 
                                    T2, ", ", Status, ") ~ ", uni.x, sep = ""))
        temp.model <- coxph(formula, data = data, method = "efron")
      }
    }
  }
  else if (!is.null(initial.model$coefficients)) {
    temp.model <- initial.model
  }
  i <- 0
  break.rule <- TRUE
  while (break.rule) {
    i <- i + 1
    if (i == 1) {
      variable.list2 <- setdiff(variable.list, all.vars(temp.model$formula))
    }
    else {
      variable.list2 <- setdiff(variable.list, c(all.vars(temp.model$formula), 
                                                 out.x))
      out.x <- NULL
    }
    if (length(variable.list2) != 0) {
      anova.pvalue <- NULL
      mv.pvalue <- NULL
      vif.value <- NULL
      for (k in 1:length(variable.list2)) {
        model <- update(temp.model, as.formula(paste(". ~ . + ", 
                                                     variable.list2[k], sep = "")))
        if (length(model$coefficients) > 1) {
          if (sum(is.na(model$coefficients)) != 0) {
            anova.pvalue[k] <- 1
            mv.pvalue[k] <- 1
            vif.value[k] <- 999
          }
          else {
            anova.pvalue[k] <- anova(temp.model, model)[2, 
                                                        "P(>|Chi|)"]
            mv.pvalue[k] <- summary(model)$coefficients[nrow(summary(model)$coefficients), 
                                                        "Pr(>|z|)"]
            model.vif <- vif(glm(as.formula(paste(Status, 
                                                  paste(names(model$coefficients), collapse = "+"), 
                                                  sep = "~")), data = data, family = binomial(link = "logit")))
            vif.value[k] <- model.vif[length(model.vif)]
          }
        }
      }
      variable.list2.1 <- variable.list2[mv.pvalue <= 0.9 & 
                                           !is.na(mv.pvalue) & vif.value <= vif.threshold]
      anova.pvalue2 <- anova.pvalue[mv.pvalue <= 0.9 & 
                                      !is.na(mv.pvalue) & vif.value <= vif.threshold]
      mv.pvalue2 <- mv.pvalue[mv.pvalue <= 0.9 & !is.na(mv.pvalue) & 
                                vif.value <= vif.threshold]
      enter.x <- variable.list2.1[anova.pvalue2 == min(anova.pvalue2, 
                                                       na.rm = TRUE) & anova.pvalue2 <= sle]
      wald.p <- mv.pvalue2[anova.pvalue2 == min(anova.pvalue2, 
                                                na.rm = TRUE) & anova.pvalue2 <= sle]
      if (length(setdiff(enter.x, NA)) != 0) {
        if (length(enter.x) > 1) {
          enter.x <- enter.x[which.min(wald.p)]
        }
        temp.model <- update(temp.model, as.formula(paste(". ~ . + ", 
                                                          enter.x, sep = "")))
      }
    }
    else {
      enter.x <- NULL
    }
    if (i == 1 & length(enter.x) == 0) {
      return(temp.model)
    }
    else {
      variable.list3 <- setdiff(rownames(summary(temp.model)$coefficients), 
                                c(enter.x, in.variable))
      if (length(variable.list3) != 0) {
        anova.pvalue <- NULL
        for (k in 1:length(variable.list3)) {
          model <- update(temp.model, as.formula(paste(". ~ . - ", 
                                                       variable.list3[k], sep = "")))
          anova.pvalue[k] <- anova(model, temp.model)[2, 
                                                      "P(>|Chi|)"]
        }
        out.x <- variable.list3[anova.pvalue == max(anova.pvalue, 
                                                    na.rm = TRUE) & anova.pvalue > sls]
        out.x <- setdiff(out.x, NA)
        if (length(out.x) != 0) {
          if (length(out.x) > 1) {
            out.x.1 <- out.x
            for (j in 1:length(out.x)) {
              out.x[j] <- out.x.1[(length(out.x) - j + 
                                     1)]
            }
            wald.p <- rep(NA, length(out.x))
            for (j in 1:length(out.x)) {
              wald.p[j] <- summary(temp.model)$coefficients[, 
                                                            "Pr(>|z|)"][rownames(summary(temp.model)$coefficients) == 
                                                                          out.x[j]]
            }
            out.x <- out.x[which.max(wald.p)]
          }
          temp.model <- update(temp.model, as.formula(paste(". ~ . - ", 
                                                            out.x, sep = "")))
        }
      }
      else {
        out.x <- NULL
      }
    }
    if ((length(enter.x) + length(out.x)) == 0) {
      final.model <- temp.model
      break.rule <- FALSE
    }
    enter.x <- NULL
  }
  return(final.model)
}

myRFLike <- function(test_exp5, 
                     test_survival4,
                     seeds = c(1,2),
                     mtry = c(10,15,20,30),
                     ntree = c(20,50,100,200)){
  re <- lapply(seeds, function(seed){
    with_seed(seed = seed, svfold_cv <- vfold_cv(test_survival4, v = 3, repeats = 5, strata = OS))
    svfold_cv <- svfold_cv[[1]]
    a <- 1
    re <- lapply(svfold_cv, function(cv){
      # Split data 
      train_survival <- training(cv)
      val_survival <- testing(cv)
      train_exp <- test_exp5[,match(train_survival$sample, colnames(test_exp5))]
      val_exp <- test_exp5[,match(val_survival$sample, colnames(test_exp5))]
      # model construction
      tuneGrid_RF <- expand.grid(mtry = mtry,
                                 ntree = ntree)
      auc1 <- lapply(1:nrow(tuneGrid_RF), function(x){
        mtry1 <- tuneGrid_RF[x,1]
        ntree1 <- tuneGrid_RF[x,2]
        # filter cox for accelerating randomForest
        cox_re <- coxFilter(train_exp, train_survival)
        cox_re2 <- cox_re[cox_re < .1]
        p_val <- .01
        while(length(cox_re2) > 25) {
          cox_re2 <- cox_re2[cox_re2 < p_val]
          p_val = p_val - 0.005
        }
        if (length(cox_re2) > 0){
          train_exp2 <- train_exp[rownames(train_exp) %in% names(cox_re2),]
        } else{
          train_exp2 <- train_exp
        }
        # Run model
        with_seed(seed = 1, rf_model <- randomForest(x=t(train_exp), y = train_survival$OS, 
                                                     importance = F, 
                                                     mtry = mtry1,
                                                     ntree = ntree1, proximity=F ))
        ## Train
        rf.pre.train <- as.data.frame(predict(rf_model, data=t(train_exp)))
        rf.pre.train$True <- train_survival$OS
        rf.pre.train <- na.omit(rf.pre.train)
        train.auc <- tryCatch(calAUC(rf.pre.train),
                              error = function(x){0})
        ## Test
        rf.pre.test <- as.data.frame(predict(rf_model, newdata=t(val_exp)))
        rf.pre.test$True <- val_survival$OS
        rf.pre.test <- na.omit(rf.pre.test)
        test.auc <- tryCatch(calAUC(rf.pre.test),
                             error = function(x){0})
        index <- a
        sel_genes <- rownames(rf_model$importance)
        sel_genes <- paste0(sel_genes, collapse = ", ")
        auc <- c(train.auc,
                 test.auc,
                 seed, index, mtry1, ntree1, sel_genes)
        # print(auc)
      })
      a <<- a + 1
      auc1 <- do.call(rbind, auc1)
      return(auc1)
    })
  })
  return(re)
}

# Validation ----
valProcess <- function(test_exp,
                       test_survival,
                       genes){
  clinical1 <- data.frame(
    patient = colnames(test_exp),
    group = ifelse(as.numeric(substr(colnames(test_exp), nchar(colnames(test_exp))-2, nchar(colnames(test_exp))-1)) < 10, "Tumor", "Normal")
  )
  clinical1$group <- factor(clinical1$group, levels = c("Normal", "Tumor"))
  # Filter
  ## DEGs
  test_exp2 <- test_exp[rownames(test_exp) %in% genes,]
  if (as.numeric(table(clinical1$group)[1]) > 4){
    DEGs <- wilcoxDEGs(test_exp2, clinical1)
    DEGs2 <- DEGs[abs(DEGs$LogFC) > 1 & DEGs$P.value < .05,]
    test_exp3 <- test_exp2[rownames(test_exp2) %in% rownames(DEGs2),]
    test_exp4 <- tumorFilter(test_exp3, test_survival)[[1]]
    test_survival4 <- tumorFilter(test_exp3, test_survival)[[2]]
  } else{
    test_exp4 <- test_exp2
    test_survival4 <- test_survival
  }
  ## Cox
  cox_re <- coxFilter(test_exp4, test_survival4)
  cox_re2 <- cox_re[cox_re < .05]
  p_val <- .01
  while(length(cox_re2) > 50) {
    cox_re2 <- cox_re2[cox_re2 < p_val]
    p_val = p_val / 2
  }
  test_exp5 <- test_exp4[rownames(test_exp4) %in% names(cox_re2),]
  ## Filtered with OS data
  test_survival5 <- na.omit(test_survival4)
  test_exp5 <- test_exp5[,match(test_survival5$sample, colnames(test_exp5))]
  # Change gene names
  rownames(test_exp5) <- gsub("-", "_", rownames(test_exp5), fixed = T)
  test_survival4 <- test_survival5
  return(list(test_survival4, test_exp5))
}

valLasso_fast <- function(tune_grids_lasso, exp, survival){
  sfExport("tune_grids_lasso", "exp", "survival")
  re <- sfSapply(1:nrow(tune_grids_lasso), function(index1){
    # index1 = 1
    seed = tune_grids_lasso[index1,1]
    cv = tune_grids_lasso[index1,2]
    alpha = tune_grids_lasso[index1,3]
    path = tune_grids_lasso[index1,5]
    index2 = which(names(pathway_list) %in% path)
    pathway = pathway_list[[index2]]
    genes = pathway@geneIds
    # Train Prepare
    test_survival4 <- valProcess(test_exp,
                                 test_survival,
                                 genes)[[1]]
    test_exp5 <- valProcess(test_exp,
                            test_survival,
                            genes)[[2]]
    # Train model
    with_seed(seed = seed, svfold_cv <- vfold_cv(test_survival4, v = 3, repeats = 5, strata = OS))
    svfold_cv <- svfold_cv[[1]]
    cv <- svfold_cv[[cv]]
    train_survival <- training(cv)
    val_survival <- testing(cv)
    train_exp <- test_exp5[,match(train_survival$sample, colnames(test_exp5))]
    val_exp <- test_exp5[,match(val_survival$sample, colnames(test_exp5))]
    # Build model
    lasso_CV <- tryCatch(with_seed(seed = 66, cv.glmnet(x=t(train_exp), 
                                                        y=train_survival$OS, nlambda = 1000,alpha = alpha)),
                         error = function(x){0})
    if(!is.list(lasso_CV)) { return(0)}
    # Validate model
    # exp <- GSE84976_exp
    exp_model <- exp[rownames(exp) %in% rownames(train_exp),]
    missed_genes <- rownames(train_exp)[!rownames(train_exp) %in% rownames(exp_model)]
    if(length(missed_genes) > 0){
      alias_df <- symbol2Alias(missed_genes, exp)
      alias_df2 <- rbind(alias_df,
                         data.frame(
                           symbol = rownames(exp_model),
                           alias = rownames(exp_model)
                         ))
      exp_model <- exp[rownames(exp) %in% c(rownames(train_exp), 
                                            alias_df$alias),]
      alias_df2 <- alias_df2[match(rownames(exp_model), alias_df2$alias),]
      rownames(exp_model) <- alias_df2$symbol
    }
    if(!all(rownames(train_exp) %in% rownames(exp_model))){
      outside.auc = 0
    } else{
      # outside_survival <- GSE84976_clinical
      outside_survival <- survival
      exp_model <- exp_model[,match(rownames(outside_survival),colnames(exp_model))]
      outside_lasso <- as.data.frame(predict(lasso_CV, newx=t(exp_model), s=lasso_CV$lambda.min))
      outside_lasso$True <- outside_survival$OS
      outside_lasso <- na.omit(outside_lasso)
      outside.auc <- tryCatch(calAUC(outside_lasso),
                              error = function(x){0})
    }
    # print(index1)
    return(outside.auc)
  })
  names(re) <- tune_grids_lasso$index
  re
}

valLasso <- function(tune_grids_lasso, exp, survival){
  re <- sapply(1:nrow(tune_grids_lasso), function(index1){
    # index1 = 1
    seed = tune_grids_lasso[index1,1]
    cv = tune_grids_lasso[index1,2]
    alpha = tune_grids_lasso[index1,3]
    path = tune_grids_lasso[index1,5]
    index2 = which(names(pathway_list) %in% path)
    pathway = pathway_list[[index2]]
    genes = pathway@geneIds
    # Train Prepare
    test_survival4 <- valProcess(test_exp,
                                 test_survival,
                                 genes)[[1]]
    test_exp5 <- valProcess(test_exp,
                            test_survival,
                            genes)[[2]]
    # Train model
    with_seed(seed = seed, svfold_cv <- vfold_cv(test_survival4, v = 3, repeats = 5, strata = OS))
    svfold_cv <- svfold_cv[[1]]
    cv <- svfold_cv[[cv]]
    train_survival <- training(cv)
    val_survival <- testing(cv)
    train_exp <- test_exp5[,match(train_survival$sample, colnames(test_exp5))]
    val_exp <- test_exp5[,match(val_survival$sample, colnames(test_exp5))]
    # Build model
    lasso_CV <- tryCatch(with_seed(seed = 66, cv.glmnet(x=t(train_exp), y=train_survival$OS, nlambda = 1000,alpha = alpha)),
                         error = function(x){0})
    if(!is.list(lasso_CV)) { return(0)}
    # Validate model
    # exp <- GSE84976_exp
    exp_model <- exp[rownames(exp) %in% rownames(train_exp),]
    missed_genes <- rownames(train_exp)[!rownames(train_exp) %in% rownames(exp_model)]
    if(length(missed_genes) > 0){
      alias_df <- symbol2Alias(missed_genes, exp)
      alias_df2 <- rbind(alias_df,
                         data.frame(
                           symbol = rownames(exp_model),
                           alias = rownames(exp_model)
                         ))
      exp_model <- exp[rownames(exp) %in% c(rownames(train_exp), 
                                            alias_df$alias),]
      alias_df2 <- alias_df2[match(rownames(exp_model), alias_df2$alias),]
      rownames(exp_model) <- alias_df2$symbol
    }
    if(!all(rownames(train_exp) %in% rownames(exp_model))){
      outside.auc = 0
    } else{
      # outside_survival <- GSE84976_clinical
      outside_survival <- survival
      exp_model <- exp_model[,match(rownames(outside_survival),colnames(exp_model))]
      outside_lasso <- as.data.frame(predict(lasso_CV, newx=t(exp_model), s=lasso_CV$lambda.min))
      outside_lasso$True <- outside_survival$OS
      outside_lasso <- na.omit(outside_lasso)
      outside.auc <- tryCatch(calAUC(outside_lasso),
                              error = function(x){0})
    }
    if(index1 %% 10 == 0) {print(index1)}
    return(outside.auc)
  })
  names(re) <- tune_grids_lasso$index
  re
}

symbol2Alias <- function(symbols, exp){
  # symbol = missed_genes
  re <- lapply(symbols, function(symbol){
    index1 <- eg2alias[eg2alias$alias_symbol %in% symbol,]$gene_id
    alias_df <- eg2alias[eg2alias$gene_id %in% index1,]
    alias <- alias_df$alias_symbol[alias_df$alias_symbol %in% rownames(exp)]
    alias <- alias[1]
    c(symbol, alias)
  })
  re <- do.call(rbind, re)
  colnames(re) <- c("symbol", "alias")
  return(as.data.frame(re))
}

# Main----
directMode <- function(test_exp,
                       test_survival,
                       genes,
                       seeds = 60:61){
  # Prepare
  clinical1 <- data.frame(
    patient = colnames(test_exp),
    group = ifelse(as.numeric(substr(colnames(test_exp), nchar(colnames(test_exp))-2, nchar(colnames(test_exp))-1)) < 10, "Tumor", "Normal")
  )
  clinical1$group <- factor(clinical1$group, levels = c("Normal", "Tumor"))
  # Filter
  ## DEGs
  test_exp2 <- test_exp[rownames(test_exp) %in% genes,]
  if (as.numeric(table(clinical1$group)[1]) > 4){
    DEGs <- wilcoxDEGs(test_exp2, clinical1)
    DEGs2 <- DEGs[abs(DEGs$LogFC) > 1 & DEGs$P.value < .05,]
    test_exp3 <- test_exp2[rownames(test_exp2) %in% rownames(DEGs2),]
    test_exp4 <- tumorFilter(test_exp3, test_survival)[[1]]
    test_survival4 <- tumorFilter(test_exp3, test_survival)[[2]]
  } else{
    test_exp4 <- test_exp2
    test_survival4 <- test_survival
  }
  ## Cox
  cox_re <- coxFilter(test_exp4, test_survival4)
  cox_re2 <- cox_re[cox_re < .05]
  p_val <- .01
  while(length(cox_re2) > 50) {
    cox_re2 <- cox_re2[cox_re2 < p_val]
    p_val = p_val / 2
  }
  test_exp5 <- test_exp4[rownames(test_exp4) %in% names(cox_re2),]
  if(nrow(test_exp5) <= 3){
    return(NA)
  } else{
    ## Filtered with OS data
    test_survival5 <- na.omit(test_survival4)
    test_exp5 <- test_exp5[,match(test_survival5$sample, colnames(test_exp5))]
    # Change gene names
    rownames(test_exp5) <- gsub("-", "_", rownames(test_exp5), fixed = T)
    # Model Part
    print(paste0("Start Train Model for: ", test_survival5$Project[1]))
    ## lasso
    time1 <- system.time(lasso.re <- myLassoLike(test_exp5, test_survival5, seeds = seeds, alphas = seq(0,1,0.1)))
    print(paste0("Lasso: ", as.numeric(time1[3])))
    lasso.re <- unlist(lasso.re, recursive = F)
    lasso.re <- do.call(rbind, lasso.re)
    lasso.re <- as.data.frame(lasso.re)
    colnames(lasso.re) <- c("Training",
                            "Testing",
                            "seed", "CV_index", "alpha", "genes")
    lasso.re$tune <- paste(lasso.re$seed,
                           lasso.re$CV_index,
                           lasso.re$alpha, sep = ":")
    lasso.re <- lasso.re[,-c(3:5)]
    lasso.re$method <- "lasso"
    ## cox
    time2 <- system.time(cox.re <- myCoxLike(test_exp5, test_survival5, seeds = seeds))
    print(paste0("Cox: ", as.numeric(time2[3])))
    cox.re <- unlist(cox.re, recursive = F)
    cox.re <- do.call(rbind, cox.re)
    colnames(cox.re) <- c("Training",
                          "Testing",
                          "seed", "CV_index", "method", "genes")
    cox.re$tune <- paste(cox.re$seed,
                         cox.re$CV_index,
                         cox.re$method, sep = ":")
    cox.re <- cox.re[,-c(3:5)]
    cox.re$method <- "cox"
    ## Random Forest
    time3 <- system.time(rf.re <- myRFLike(test_exp5, test_survival5, seeds = seeds))
    print(paste0("Random Forest: ", as.numeric(time3[3])))
    rf.re <- unlist(rf.re, recursive = F)
    rf.re <- do.call(rbind, rf.re)
    rf.re <- as.data.frame(rf.re)
    colnames(rf.re) <- c("Training",
                         "Testing",
                         "seed", "CV_index", "mtry",
                         "ntree", "genes")
    rf.re$tune <- paste(rf.re$seed,
                        rf.re$CV_index,
                        rf.re$mtry,
                        rf.re$ntree, sep = ":")
    rf.re <- rf.re[,-c(3:6)]
    rf.re$method <- "randomforest"
    # merge
    model.re <- rbind(rf.re,
                      cox.re,
                      lasso.re)
    model.re$Type <- test_survival5$Project[1]
    print("Finish!")
    return(model.re)
  }
}


