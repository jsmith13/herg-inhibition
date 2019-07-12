## required libraries
require(ggplot2)
require(dplyr)
require(pROC)
require(shiny)
require(kableExtra)


## define Shiny UI
fluidPage(
  
  # a title for the page
  titlePanel("hERG Inhibition: Model Exploration"),
  
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
      
      # check boxes for higher order terms
      checkboxGroupInput("higher.order.predictors",
                         h4("Higher Order Terms"),                 
                         choices = list("VABC Volume^2" = "I(VABC.Volume.Descriptor^2)",
                                        "TPSA^2" = "I(Topological.Polar.Surface.Area^2)",
                                        "XLogP^3" = "I(XLogP^3)"),
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
        column(6,
               # display the prediction plot
               plotOutput("prediction.plot", height = "450px")
        ),
        column(4,
               # display the ROC plot
               plotOutput("roc.plot", height = "450px")
        )
      ),
      
      fluidRow(
        column(6,
               # display the modeling dataset plot
               plotOutput("modeling.dataset.plot", height = "450px")
        ),
        
        column(2,
               # display the confusion matrix
               tableOutput("confusion.matrix")
        ),
        column(2,
               # explicitly print the sensitivity and specificity
               tableOutput("sensitivity_specificity")
        )
      )
    )
  )
)