### Shiny application for visual exploration of model fitting.

## setup dataset to be used in app
# set random seed for reproducability
set.seed(18438)

# import required libraries
require(ggplot2)
require(dplyr)
require(tidyr)
require(robustbase)
require(pROC)

# import a 500 row sample from the hERG dataset
hERG <- read.csv("https://gist.githubusercontent.com/jsmith13/dd405d28a3dfa9e2a7ef276cc6135cbc/raw/23f20c56c44bd9cd7f212a745df8295ca177f6a3/hERG.csv")

# generate principle components from the major numeric descriptors
# log transform predictors where necessary to get closer to normal
hERG.PCA <- prcomp(formula = ~ log(Bond.Polarizabilities) + log(VABC.Volume.Descriptor) + log(Topological.Polar.Surface.Area + 1) + log(Molecular.Weight) + XLogP + log(Eccentric.Connectivity.Index), data = hERG, center = TRUE, scale. = TRUE)

# ten-fold partition the hERG dataset
logit.folds <- caret::createFolds(hERG$Phenotype, k = 10)



## Shiny user interface elements
UI <- shinyUI(fluidPage(
  
  # a title for the page
  titlePanel("A Logistic Regression Model for hERG Inhibition"),
  
  # lay out the page with a minor and major panel
  sidebarLayout(
    
    # a minor panel for selecting modeling options
    sidebarPanel(
      
      h3("Model Parameters"),
      
      # check boxes for the available predictors
      checkboxGroupInput("selected.predictors",
        h4("Predictors to Include"),
        choices = list("Aromatic Atoms" = "Aromatic.Atoms.Count",
                       "Bond Polarizabilities" = "Bond.Polarizabilities",
                       "VABC Volume" = "VABC.Volume.Descriptor",
                       "Hydrogen Bond Acceptors" = "Hydrogen.Bond.Acceptors",
                       "Hydrogen Bond Donors" = "Hydrogen.Bond.Donors",
                       "Rotatable Bonds" = "Rotatable.Bonds.Count",
                       "TPSA" = "Topological.Polar.Surface.Area",
                       "Molecular Weight" = "Molecular.Weight",
                       "XLogP" = "XLogP",
                       "Formal Charge" = "Formal.Charge",
                       "sp3 Character" = "SP3.Character",
                       "Eccentric Connectivity Index" = "Eccentric.Connectivity.Index",
                       "Amine Present" = "Amine"),
        selected = c("Aromatic.Atoms.Count", "VABC.Volume.Descriptor", "Amine", "Topological.Polar.Surface.Area", "XLogP")
        ),
      
      # check boxes for quadratic terms
      checkboxGroupInput("higher.order.predictors",
        h4("Higher Order Terms"),                 
        choices = list("VABC Volume" = "I(VABC.Volume.Descriptor^2)",
                       "TPSA" = "I(Topological.Polar.Surface.Area^2)",
                       "XLogP" = "I(XLogP^3)"),
        selected = c("I(VABC.Volume.Descriptor^2)", "I(Topological.Polar.Surface.Area^2)", "I(XLogP^3)")
      ),
      
      h3("Plotting Options"),
      
      # radio buttons to change plotted color scheme from binary active/inactive to odds active
      radioButtons("plot.color.scheme", h4("Prediction Plot"), choices = c("Binary Active/Inactive" = "binary", "Odds Active" = "continuous", "Correct/Incorrect" = "correct")),
      
      # a slider to select the odds cutoff for classifying active/inactive predictions
      sliderInput("threshold", h4("Active/Inactive Cutoff"), min = 0.1, max = 0.9, value = 0.5, step = 0.05),
      
      width = 2
    ),
    
    # the major panel for displaying plots and summaries
    mainPanel(
      fluidRow(
        column(12,
               # display the prediction plot
               plotOutput("prediction.plot", height = "600px")
        )
      ),
      
      fluidRow(
        column(8,
          # display the modeling dataset plot
          plotOutput("modeling.dataset.plot")
        ),
        
        column(4,
          # display the ROC plot
          plotOutput("roc.plot")     
        )
      )
      
    )
  )
))



## Shiny server
server <- shinyServer(function(input, output) {
  
  # define the modeling dataset plot
  # plot the chemical space using the principle components, color by activity score
  # set axis labels and legend parameters
  output$modeling.dataset.plot <- renderPlot(ggplot(data = hERG) + geom_point(aes(x = hERG.PCA$x[,1], y = hERG.PCA$x[,2], color = as.factor(Phenotype))) +
    ggtitle("Modeling Dataset") + xlab("Chemical Space Projection") + ylab("") + theme(axis.ticks = element_blank(), axis.text = element_blank()) + 
    guides(colour = guide_legend(title = "", reverse = TRUE)) + scale_color_manual(values = c("dark blue", "red"), labels = c("Inactive", "Active")))
  
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
    prediction.plot <- ggplot(data = hERG) + geom_point(aes(x = hERG.PCA$x[,1], y = hERG.PCA$x[,2], color = binary.continuous.switch())) +
    ggtitle("Predicted Values") + xlab("Chemical Space Projection") + ylab("") + theme(axis.ticks = element_blank(), axis.text = element_blank()) + 
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
  })
  
  # define the ROC plot
  output$roc.plot <- renderPlot({
    # calculate the ROC object
    roc.calculation <- roc(response = hERG$Phenotype, pred = predictions())

    # generate a plot; add a title and caption with the AUC
    ggroc(roc.calculation) + ggtitle("10-fold Cross-Validation ROC Curve") + labs(caption = paste("AUC:", round(roc.calculation$auc, digits = 2)))
  })
  
})



# run the Shiny app
shinyApp(ui = UI, server = server)
