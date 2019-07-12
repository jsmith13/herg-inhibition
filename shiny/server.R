## required libraries
require(ggplot2)
require(dplyr)
require(pROC)
require(shiny)
require(kableExtra)


## setup dataset to be used in app
# import a 500 row sample from the hERG dataset
hERG <- read.csv("https://gist.githubusercontent.com/jsmith13/dd405d28a3dfa9e2a7ef276cc6135cbc/raw/23f20c56c44bd9cd7f212a745df8295ca177f6a3/hERG.csv")

# generate principle components from the major numeric descriptors
# log transform predictors where necessary to get closer to normal
hERG.PCA <- prcomp(formula = ~ log(Bond.Polarizabilities) + log(VABC.Volume.Descriptor) + log(Topological.Polar.Surface.Area + 1) + log(Molecular.Weight) + XLogP + log(Eccentric.Connectivity.Index), data = hERG, center = TRUE, scale. = TRUE)

# set random seed for reproducability
set.seed(18438)

# ten-fold partition the hERG dataset
logit.folds <- caret::createFolds(hERG$Phenotype, k = 10)

## define Shiny server
function(input, output, session) {

  # render the plot of observed values
  output$observed.plot <- renderPlot(
    {
      # define the plot
      # plot the chemical space using the principle components, color by activity score
      # set axis labels and legend parameters  
      ggplot(data = hERG) + 
      geom_point(aes(x = hERG.PCA$x[,1], y = hERG.PCA$x[,2], color = as.factor(Phenotype))) +
      ggtitle("Observed Classes") + xlab("Chemical Space Projection") + ylab("") + 
      theme(axis.ticks = element_blank(), axis.text = element_blank()) + 
      guides(colour = guide_legend(title = "", reverse = TRUE)) + 
      scale_color_manual(values = c("dark blue", "red"), labels = c("Inactive", "Active"))
    },
    
    # define a function to specify the height of the plot from the width of the session
    height = function() {
      0.75 * session$clientData$output_observed.plot_width
    }
  )
  
  # define a reactive expression that reads the selected.predictors check box group and outputs a function
  selected.formula <- eventReactive(c(input$selected.predictors, input$higher.order.predictors), {
    # paste the selected terms into a formula
    # add the higher order terms if any are selected
    if (length(input$higher.order.predictors) > 0) {
      selected.formula <- paste("Phenotype ~", paste(input$selected.predictors, collapse = " + "), " + ", paste(input$higher.order.predictors, collapse = " + ")) %>% as.formula()
    }
    # exclude higher order terms otherwise
    else {
      selected.formula <- paste("Phenotype ~", paste(input$selected.predictors, collapse = " + "), collapse = " + ") %>% as.formula()
    }
    
    # return the formula
    return(selected.formula)
  })
  
  # define a reactive expression outputs predictions from a 10-fold cross-validation
  predictions <- reactive({
    # declare a dummy vector to hold the predicted values
    predictions <- 1:length(hERG$Phenotype)
    
    # perform a ten-fold cross-validation to generate predicted values 
    for (i in 1:10) {
      # fit the models on the 9/10 sets
      glm.part <- glm(selected.formula(), data = hERG[-logit.folds[[i]], ], family = binomial)
      
      # predict with the 1/10 sets
      pred.part <- predict(glm.part, hERG[logit.folds[[i]], ])
      
      # store predicted values
      predictions[logit.folds[[i]]] <- pred.part
    }
    
    # return predictions
    return(predictions)
  })
  
  # define a reactive expression that reads the plot.color.scheme radio buttons and outputs the appropriate vector
  binary.continuous.switch <- reactive({
    if (input$plot.color.scheme == "binary"){
      # calculate the odds active and cut into two groups at the user defined threshold
      return(boot::inv.logit(predictions()) %>% cut(breaks = c(-0.01, input$threshold, 1)))
    }
    else if (input$plot.color.scheme == "continuous") {
      # calculate the odds active
      return(boot::inv.logit(predictions()))
    }
    else if (input$plot.color.scheme == "correct") {
      # generate a vector of active/inactive predictions
      predicted.phenotype <- boot::inv.logit(predictions()) %>% cut(breaks = c(-0.01, input$threshold, 1))
      
      # compare to known phenotype to generate a vector classifying predictions as correct/false
      return(hERG$Phenotype == as.numeric(predicted.phenotype) - 1)
    }
  })
  
  # define the prediction plot
  output$prediction.plot <- renderPlot({
    # define the cross-validation prediction plot
    # plot the chemical space using the principle components, color by predicted odds active
    # set axis labels 
    prediction.plot <- ggplot(data = hERG) + 
      geom_point(aes(x = hERG.PCA$x[,1], y = hERG.PCA$x[,2], color = binary.continuous.switch())) +
      ggtitle("Predicted Classes") + xlab("Chemical Space Projection") + ylab("") + 
      theme(axis.ticks = element_blank(), axis.text = element_blank()) + 
      guides(colour = guide_legend(title = "", reverse = TRUE))
    
    # set the legend and color scheme for the continuous and binary output options
    if (input$plot.color.scheme == "binary") {
      prediction.plot + scale_color_manual(values = c("dark blue", "red"), labels = c("Inactive", "Active"))
    }
    else if (input$plot.color.scheme == "continuous") {
      prediction.plot + scale_color_gradient("Odds Active", low = "dark blue", high = "red")
    }
    else if (input$plot.color.scheme == "correct") {
      prediction.plot + scale_color_manual(values = c("orange2", "dark blue"), labels = c("Incorrect", "Correct"))
    }
  },
    
    # define a function to specify the height of the plot from the width of the session
    height = function() {
      0.75 * session$clientData$output_prediction.plot_width
    }
  )
  
  # define the confusion matrix
  output$confusion.matrix <- renderPrint({
    # convert predicted log(Pr) into probability active and cut at user-defined threshold
    predicted.classes <- boot::inv.logit(predictions()) %>% cut(breaks = c(-0.01, input$threshold, 1)) %>% as.numeric()
    
    # generate a confusion matrix
    cm <- table(observed = hERG$Phenotype, predicted = predicted.classes - 1)
    
    # give row and column names to the matrix
    cm <- as.matrix(cm)
    rownames(cm) <- c("Inactive", "Active")
    colnames(cm) <- c("Inactive", "Active")
    
    # format the table for display
    knitr::kable(cm, caption = "Confusion Matrix") %>% kable_styling(c("bordered", "hover"), full_width = FALSE)
  })
  
  # define the sensitivity and specificity table
  output$sensitivity_specificity <- renderPrint({
    # convert predicted log(Pr) into probability active and cut at user-defined threshold
    predicted.classes <- boot::inv.logit(predictions()) %>% cut(breaks = c(-0.01, input$threshold, 1)) %>% as.numeric()
    
    # generate a confusion matrix
    cm <- table(observed = hERG$Phenotype, predicted = predicted.classes - 1)
    
    # calculate the sensitivity, specificity, and positive predictive value
    sensitivity <- round(cm[4] / (cm[2] + cm[4]), 2)
    specificity <- round(cm[1] / (cm[1] + cm[3]), 2)
    ppv <- round(cm[4] / (cm[3] + cm[4]), 2)
    
    #calculate the ROC AUC
    roc.calculation <- roc(response = hERG$Phenotype, pred = predictions())
    auc <- round(roc.calculation$auc, 2)
    
    # put into a table to print
    c("Sensitivity", "Specificity", "PPV", "AUC", sensitivity, specificity, ppv, auc) %>%
      matrix(ncol = 2) %>% kable(caption = "Performance") %>% kable_styling(c("bordered", "hover"), full_width = FALSE)
  })
  
  # define the ROC plot
  output$roc.plot <- renderPlot(
    {
      ## for plotting the ROC curve
      # calculate the ROC object
      roc.calculation <- roc(response = hERG$Phenotype, pred = predictions())
      
      # generate a plot; add a title and caption with the AUC
      to_plot <- ggroc(roc.calculation, size = 1.5) + ggtitle("10-fold Cross-Validation ROC Curve") + theme_bw()
      
      ## for highlighting the effect of the current threshold
      # convert predicted log(Pr) into probability active and cut at user-defined threshold
      predicted.classes <- boot::inv.logit(predictions()) %>% cut(breaks = c(-0.01, input$threshold, 1)) %>% as.numeric()
      
      # generate a confusion matrix
      cm <- table(observed = hERG$Phenotype, predicted = predicted.classes - 1)
      
      # calculate the sensitivity and specificity
      sens <- cm[4] / (cm[2] + cm[4])
      spec <- cm[1] / (cm[1] + cm[3])
      
      # add a vertical line at the current specificity to the ROC curve
      to_plot + geom_segment(aes(x = spec, y = sens, xend = spec, yend = 0), size = 1, color = "blue") +
        geom_segment(aes(x = spec, y = sens, xend = 1, yend = sens), size = 1, color = "blue")
    },
    
    # define a function to specify the height of the plot from the width of the session
    height = function() {
      session$clientData$output_roc.plot_width
    }
  )
  
}