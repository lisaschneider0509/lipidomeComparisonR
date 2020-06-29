#####PIMA INDIANS DATASET

##load libraries-------------------------------------------------------------------------------------------------
library(shiny)
library(shinythemes)
library(utils)
library(gridExtra)
library("MASS")
library(ggplot2) 
library(car)

pima <- rbind(MASS::Pima.te, MASS::Pima.tr)

##definition for scatterplot-------------------------------------------------------------------------------------
panel.cor <- function(x,y,digits = 2,prefix = "",cex.cor,...)
{ usr <- par("usr")
on.exit(par(usr))
par(usr = c(0, 1, 0, 1))
r <- abs(cor(x, y))
txt <- format(c(r, 0.123456789), digits = digits)[1]
txt <- paste0(prefix, txt)
if (missing(cex.cor))
  cex.cor <- 0.8 / strwidth(txt)
text(0.5, 0.5, txt, cex = cex.cor * 0.5)
}

#----------------------------------------------------###UI###----------------------------------------------------#

ui <- fluidPage(theme = shinytheme("flatly"),
                navbarPage(
                  tags$h2(p(strong('Pima'))),
                  
                  ##tab panel - data overview--------------------------------------------------------------------------------------  
                  tabPanel(tags$h3(p(strong('Data Overview'))), 
                           tags$h2("Data Overview"), 
                           hr(), 
                           mainPanel( 
                             tabsetPanel(
                               tabPanel("Description", verbatimTextOutput("description")),
                               tabPanel("Rawdata", tableOutput("rawdataPima"))
                             ))),
                  
                  #tab panel - data exploration---------------------------------------------------------------------------------
                  tabPanel(tags$h3(p(strong('Data Exploration'))),
                           tags$h3("Distribution of Pima Indians datasets "), hr(), br(),
                           
                           sidebarLayout(sidebarPanel(radioButtons(
                             "dataset",
                             "Select a Dataset", choices = c(
                               "Nr of pregn",
                               "glucose",
                               "blood pressure",
                               "skin",
                               "bmi",
                               "ped",
                               "age", 
                               "Type [Yes/No]"))),
                             
                             mainPanel(tags$h4("Data and Visualization"),
                                       tabsetPanel(
                                         tabPanel("Summary", verbatimTextOutput("summary")),
                                         tabPanel("QQ-Plot", plotOutput("qqplot")),
                                         tabPanel("Histogram", plotOutput("hist")),
                                         tabPanel("Boxplot", plotOutput("boxplot")),
                                         tabPanel("Barplot", plotOutput("barplot"))
                                       )))),
                  
                  ##tab panel - Correlation-------------------------------------------------------------------------------------
                  tabPanel(tags$h3(p(strong('Correlation'))), 
                           tags$h2("Overview Correlation"),
                           hr(),
                           mainPanel(plotOutput("scatter"))
                  ),
                  
                  ##tab panel - logistic model------------------------------------------------------------------------------------
                  tabPanel(tags$h3(p(strong('Regression Model'))), 
                           tags$h3("Choose model settings"),
                           hr(),
                           sidebarLayout(
                             sidebarPanel(
                               selectInput("regressand", "Dependent Variable", choices = c("npreg", "glu", "bp", "skin", "bmi", "ped", "age", "type")),
                               checkboxGroupInput("independent_var", 
                                                  "Independent variables",
                                                  choiceNames = c("Nr of pregnancies", "plasma glucose conc", 
                                                                  "blood pressure", "skin fold thickness", "BMI", "ped", "age"),
                                                  choiceValues = c("npreg", "glu", 
                                                                   "bp", "skin", "bmi", "ped", "age"),
                                                  selected = c("npreg", "glu", 
                                                               "bp", "skin", "bmi", "ped", "age")),
                               radioButtons("transformation", 
                                            "Transformations", 
                                            choices = c("Without Transformation", 
                                                        "Log(X)", 
                                                        "Log(Y)", 
                                                        "Log/Log", 
                                                        "X^2", 
                                                        "Y^2", 
                                                        "X^2/Y^2")), 
                               checkboxGroupInput("checkGroup", 
                                                  label = h4("Remove Outliers"), 
                                                  choices = c(rownames(pima)),  
                                                  selected = c(rownames(pima)))),
                             mainPanel(tags$h4("Regression Models"), hr(),
                                       tabsetPanel(
                                         tabPanel("Step (AIC)", 
                                                  tags$h4("AIC"), verbatimTextOutput("stepmodel")),
                                         tabPanel("Residue Plots", 
                                                  tags$h4("Theorem"), verbatimTextOutput("residue"),
                                                  tags$h4("Residue Plots"), plotOutput("residuePlot")),
                                         tabPanel("Model Summary", 
                                                  tags$h4("Theorem"), verbatimTextOutput("modelTheorem"),
                                                  tags$h4("Model Summary"), verbatimTextOutput("modelSummary"))
                                         
                                         # ,tabPanel("Final Model & Plot", 
                                         #          #tags$h4("Model equation"), verbatimTextOutput("finalModel"), 
                                         #          tags$h4("Regression curve", plotOutput("modelPlot")))
                                       ))
                           ))
                ))

### SERVER ### ----------------------------------------------------------------------------------------
server <- function(input, output) {
  ## datasets ##
  datasetInput <- reactive({
    switch(input$dataset,
           "Nr of pregn" = pima[1], 
           "glucose" = pima[2],
           "blood pressure" = pima[3],
           "skin" = pima[4],
           "bmi" = pima[5],
           "ped" = pima[6],
           "age" = pima[7],
           "Type [Yes/No]" = pima[8])})
  
  datasetNames <- reactive({
    switch (input$dataset,
            "Nr of pregn" = "Number of Pregnancies", 
            "glucose" = "Glucose",
            "blood pressure" = "Blood pressure (mm Hg)",
            "skin" = "Skin (mm)",
            "bmi" = "BMI",
            "ped" = "Diabetes pedigree function",
            "age" = "Age",
            "Type [Yes/No]" = "Type [Yes/No]")})
  
  datasetRegressand <- reactive({
    switch (input$regressand,
            "npreg" = "npreg", 
            "glu" = "glu", 
            "bp"= "bp", 
            "skin" = "skin", 
            "bmi" = "bmi", 
            "ped" = "ped", 
            "age" = "age", 
            "type" = "type")})
  
  ## dataoverview -------------------------------------------------------------------------------------------------------
  output$rawdataPima <- renderTable({
    # print rawdata
    dataset <- pima
    dataset})
  
  output$description <- renderText({  #print hepppage
    "Pima.tr {MASS}	    R Documentation
    
    Diabetes in Pima Indian Women
    
    Description:
    A population of women who were at least 21 years old, of Pima Indian heritage and living near Phoenix, Arizona, was tested for diabetes according to World Health Organization criteria. 
    The data were collected by the US National Institute of Diabetes and Digestive and Kidney Diseases. We used the 532 complete records after dropping the (mainly missing) data on serum insulin.
    
    Usage:
    Pima.tr
    Pima.tr2
    Pima.te
    
    Format:
    npreg ...     number of pregnancies.
    glu ...       plasma glucose concentration in an oral glucose tolerance test.
    bp ...        diastolic blood pressure (mm Hg).
    skin ...      triceps skin fold thickness (mm).
    bmi ...       body mass index (weight in kg/(height in m)^2).
    ped ...       diabetes pedigree function.
    age ...       age in years.
    type ...      Yes or No, for diabetic according to WHO criteria."
  })
  
  ## Data Exploration -------------------------------------------------------------------------------------------------------
  output$summary <- renderPrint({
    dataset <- datasetInput()
    summary(dataset[,1])})
  
  output$qqplot <- renderPlot({
    dataset <- datasetInput()
    caption <- datasetNames()
    if(caption == "Type [Yes/No]"){}
    else{
      qqnorm(dataset[,1], main = caption)
      qqline(dataset, col=2)}})
  
  output$hist <- renderPlot({
    dataset <- datasetInput()
    caption <- datasetNames()
    if(caption == "Type [Yes/No]"){}
    else{
      hist(dataset[,1], freq = FALSE, main = caption, xlab = caption)
      lines(density(dataset[,1]), lwd=2)
      lines(density(dataset[,1],adjust=1.5),col=2,lwd=2)}
    # ggplot(pima, aes(x=dataset[,1])) + geom_histogram(binwidth = 1, 
    # aes(y= ..density.., fill = ..count..))+geom_density(fill="coral", alpha = 0.4)   + labs(x="")
  })
  
  output$boxplot <- renderPlot({
    dataset <- datasetInput()
    caption <- datasetNames()
    if(caption == "Type [Yes/No]"){}
    else{
      Boxplot(y = dataset[, 1], main = caption, ylab = caption)}})
  
  output$barplot <- renderPlot({
    dataset <- datasetInput()
    caption <- datasetNames()
    if(caption == "Type [Yes/No]"){
      barplot(table(dataset[,1]), main = caption)}
    else{}
  })
  
  ## correlation -------------------------------------------------------------------------------------------------------
  output$scatter <- renderPlot({
    dataset <- datasetInput()
    pairs(pima, lower.panel = panel.smooth, upper.panel = panel.cor,
          gap=0, row1attop=FALSE, main = "Scatterplot")})
  
  
  ## logistisches model--------------------------------------------------------------------------------------------------  
  regModel <- reactive({
    if (input$transformation == "Without Transformation") 
    {expln <- paste(input$independent_var, collapse = "+")
    as.formula(paste(input$regressand, "~", expln))} 
    else if (input$transformation == "Log(X)") 
    {expln <- paste("log(", input$independent_var, ")", collapse = "+")
    as.formula(paste(input$regressand, "~", expln))} 
    else if (input$transformation == "Log(Y)")
    {expln <- paste(input$independent_var, collapse = "+")
    as.formula(paste("log(",input$regressand, ")", "~", expln))}
    else if (input$transformation == "Log/Log")
    {expln <- paste("log(", input$independent_var, ")", collapse = "+")
    as.formula(paste("log(",input$regressand, ")", "~", expln))}
    else if (input$transformation == "Standardisation")
    {expln <- paste(input$independent_var, collapse = "+")
    as.formula(paste(input$regressand, "~", expln))} 
    else if (input$transformation == "X^2") {
      {expln <- paste(input$independent_var, "^2", collapse = "+")
      as.formula(paste(input$regressand, "~", expln))}
    }
    else if (input$transformation == "Y^2") {
      {expln <- paste(input$independent_var, collapse = "+")
      as.formula(paste(input$regressand, "^2", "~", expln))}
    }
    else if (input$transformation == "X^2/Y^2"){
      expln <- paste(input$independent_var, "^2", collapse = "+")
      as.formula(paste(input$regressand, "^2", "~", expln))}
  })
  
  ## @AIC ------------------------------
  output$AIC  <- renderPrint({
    linearModel()})
  
  output$stepmodel <- renderPrint({
    reg <- datasetRegressand()
    if(reg == "type"){
      fit = glm(regModel(), family = binomial(link = "logit"), data=pima[c(input$checkGroup),])
      step(fit)
    }
    else{
      fit = lm(regModel(), data=pima[c(input$checkGroup),])
      step(fit)}
  })
  
  ## residue plots
  output$residue <- renderPrint({
    regModel()})
  
  output$residuePlot <- renderPlot ({
    reg <- datasetRegressand()
    if(reg == "type"){
      fit = glm(regModel(), family = binomial(link = "logit"), data=pima[c(input$checkGroup),])
      par(mfrow=c(2,2))
      plot(fit)
    }
    else{
      fit = lm(regModel(), data=pima[c(input$checkGroup),])
      par(mfrow=c(2,2))
      plot(fit)}
  })
  
  ## Model Summary
  output$modelTheorem <- renderPrint({
    regModel()})
  
  output$modelSummary <- renderPrint({
    reg <- datasetRegressand()
    if(reg == "type"){
      fit = glm(regModel(), family = binomial(link = "logit"), data=pima[c(input$checkGroup),])
      par(mfrow=c(2,2))
      summary(fit)
    }
    else{
      fit = lm(regModel(), data=pima[c(input$checkGroup),])
      par(mfrow=c(2,2))
      summary(fit)}
  })
  
  # output$finalModel <- renderText({
  #   fit = glm(regModel(), family = binomial(link = "logit"), data=pima[c(input$checkGroup),])
  #   fit <- summary(fit)
  #   fit <- fit$coefficients
  #   alpha <- fit[1]
  #   beta <- fit[2]
  #   sprintf("y = %s + %s * x", alpha, beta)})
  
  # output$modelPlot <- renderPlot({
  #   fit = glm(regModel(), family = binomial(link = "logit"), data=pima[c(input$checkGroup),])
  #   {plot(x = c(1, 2, 3), y = c(6, 5, 4))
  #   abline(fit, col = 2, lwd = 3)}
  # })
  #   
  }

shinyApp(ui = ui, server = server)