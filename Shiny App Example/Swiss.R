#####SWISS DATASET

##load libraries-------------------------------------------------------------------------------------------------
library(shiny)
library(shinythemes)
library(utils)
library(gridExtra)
library(ggplot2) 
library(plyr)
library(car)

##normalize data before analysis---------------------------------------------------------------------------------
Fertility <- scale(swiss$Fertility, center = TRUE, scale = TRUE)
Agriculture <- scale(swiss$Agriculture, center = TRUE, scale = TRUE)
Examination <- scale(swiss$Examination, center = TRUE, scale = TRUE)
Education <- scale(swiss$Education, center = TRUE, scale = TRUE)
Infant.Mortality <- scale(swiss$Infant.Mortality, center = TRUE, scale = TRUE)
Catholic <- scale(swiss$Catholic, center = TRUE, scale = TRUE)
normSwiss <- cbind(Fertility, Agriculture, Examination, Education, Catholic, Infant.Mortality)
colnames(normSwiss) <- c("Fertility", "Agriculture", "Examination", "Education", "Catholic", "Infant.Mortality")
rownames(normSwiss) <- c(row.names(swiss))
normSwiss <- as.data.frame(normSwiss)

##definition for scatterplot-------------------------------------------------------------------------------------
panel.cor <- function(x,y,digits = 2,prefix = "",cex.cor,...)
{
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if (missing(cex.cor))
    cex.cor <- 0.8 / strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * 0.5)
}

#-------------------------------------------------###UI###-------------------------------------------------------#

ui <- fluidPage(theme = shinytheme("flatly"),
                
                navbarPage(
                  tags$h2(p(strong('Swiss Data'))),
                  
                  ##tab panel - data overview----------------------------------------------------------------------------------- 
                  tabPanel(tags$h3(p(strong('Data Overview'))), 
                           tags$h2("Data Overview"), 
                           hr(), 
                           mainPanel( 
                             tabsetPanel(
                               tabPanel("Description", verbatimTextOutput("description")),
                               tabPanel("Rawdata", tableOutput("rawdataSwiss"))
                             ))),
                  
                  ##tab panel - data exploration--------------------------------------------------------------------------------
                  tabPanel(tags$h3(p(strong('Data Exploration'))),
                           tags$h3("Data Exploration: Distribution of Swiss Dataset "), hr(), br(),
                           
                           sidebarLayout(sidebarPanel(radioButtons(
                             "dataset",
                             "Please select a dataset",choices = c(
                               "Fertility",
                               "Agriculture",
                               "Education",
                               "Catholic",
                               "Infant.Mortality"))),
                             
                             mainPanel(tags$h4("Data and Visualization"),
                                       tabsetPanel(
                                         tabPanel("Summary", 
                                                  h4(textOutput("caption")),
                                                  verbatimTextOutput("summary")),
                                         tabPanel("QQ-Plot", 
                                                  plotOutput("qqplot")),
                                         tabPanel("Histogram", 
                                                  plotOutput("hist")),
                                         tabPanel("Boxplot", 
                                                  plotOutput("boxplot")
                                         )
                                       )))),
                  
                  ##tab panel - correlation-------------------------------------------------------------------------------------
                  tabPanel(tags$h3(p(strong('Correlation'))), 
                           tags$h2("Overview Correlation"),
                           hr(),
                           mainPanel(plotOutput("scatter"))),
                  
                  ##tab panel - linear model------------------------------------------------------------------------------------
                  tabPanel(tags$h3(p(strong('Linear Model'))), 
                           tags$h3("Please choose dependent and independent Variables"),
                           hr(),
                           sidebarLayout(
                             sidebarPanel(
                               selectInput("regressand", 
                                           "Dependent Variable", 
                                           choices = c("Fertility", 
                                                       "Agriculture", 
                                                       "Examination", 
                                                       "Education", 
                                                       "Catholic", 
                                                       "Infant.Mortality" )),
                               checkboxGroupInput("checkbox", 
                                                  label = h5(strong("Independent variables")), 
                                                  choices = c(colnames(swiss)), 
                                                  # choiceValues = c("Fertility", 
                                                  #                  "Agriculture", 
                                                  #                  "Examination", 
                                                  #                  "Education", 
                                                  #                  "Catholic", 
                                                  #                  "Infant.Mortality"), 
                                                  selected = c(colnames(swiss))),
                               radioButtons("transformation", 
                                            "Apply", 
                                            choices = c("Without Transformation", 
                                                        "Standardisation", 
                                                        "Log(X)", 
                                                        "Log(Y)", 
                                                        "Log/Log", 
                                                        "X^2", 
                                                        "Y^2", 
                                                        "X^2/Y^2")),
                               checkboxGroupInput("checkGroup", 
                                                  label = h5(strong("Select data/remove outliers")), 
                                                  choices = c(rownames(swiss)),  
                                                  selected = c(rownames(swiss)))
                             ),
                             mainPanel(tags$h4("Linear Models"), hr(),
                                       tabsetPanel(
                                         tabPanel("Step (AIC)", 
                                                  tags$h4("AIC"), verbatimTextOutput("stepmodel")),
                                         tabPanel("Residue Plots", 
                                                  tags$h4("Theorem"), verbatimTextOutput("modelEquation"),
                                                  tags$h4("Residue Plots"), plotOutput("residuePlot")),
                                         tabPanel("Model Summary", 
                                                  tags$h4("Theorem"), verbatimTextOutput("modelTheorem"),
                                                  tags$h4("Model Summary"), verbatimTextOutput("modelSummary"))
                                         # , tabPanel("Final Model & Plot",
                                         #          #tags$h4("Model equation"), verbatimTextOutput("finalModel"),
                                         #          tags$h4("Regression curve"), verbatimTextOutput("modelPlot"))
                                       ))
                           ))))

#--------------------------------------------------###SERVER###--------------------------------------------------#

server <- function(input, output){
  datasetInput <- reactive({
    switch(input$dataset,
           "Fertility" = swiss[1], 
           "Agriculture" = swiss[2],
           "Examination" = swiss[3],
           "Education" = swiss[4],
           "Catholic" = swiss[5],
           "Infant.Mortality" = swiss[6])})
  
  datasetNames <- reactive({
    switch(input$dataset,
           "Fertility" = "Fertility", 
           "Agriculture" = "Agriculture",
           "Examination" = "Examination",
           "Education" = "Education",
           "Catholic" = "Catholic",
           "Infant.Mortality" = "Infant.Mortality")})
  
  
  ##description-------------------------------------------------------------------------------------------------- 
  output$caption<-renderText({
    switch(input$dataset,
           "Fertility" = "Fertility", 
           "Agriculture" = "Agriculture",
           "Examination" = "Examination",
           "Education" = "Education",
           "Catholic" = "Catholic",
           "Fertility" = "Fertility",
           "Infant.Mortality" = "Infant.Mortality")})
  
  output$value <- renderText({output$caption})
  
  ##rawdata--------------------------------------------------------------------------------------------------------
  output$rawdataSwiss <- renderTable({
    dataset <- swiss
    dataset})  #head(dataset)
  
  output$description <- renderText({
    "swiss {datasets}	                    R Documentation
    
    Swiss Fertility and Socioeconomic Indicators (1888) Data
    
    Description:
    Standardized fertility measure and socio-economic indicators for each of 47 French-speaking provinces of Switzerland at about 1888.
    
    Usage:
    swiss
    
    Format: 
    A data frame with 47 observations on 6 variables, each of which is in percent, i.e., in [0, 100].
    [,1] 	Fertility 	Ig, ‘common standardized fertility measure’
    [,2] 	Agriculture	% of males involved in agriculture as occupation
    [,3] 	Examination	% draftees receiving highest mark on army examination
    [,4] 	Education 	% education beyond primary school for draftees.
    [,5] 	Catholic 	% ‘catholic’ (as opposed to ‘protestant’).
    [,6] 	Infant.Mortality	live births who live less than 1 year. "
  })
  
  ##summary--------------------------------------------------------------------------------------------------------
  output$summary <- renderPrint({
    dataset <- datasetInput()
    summary(dataset[, 1])})
  
  ##histogram--------------------------------------------------------------------------------------------------------
  output$hist <- renderPlot({
    dataset <- datasetInput()
    header <- datasetNames()
    ggplot(swiss, aes(x=dataset[, 1])) + geom_histogram(binwidth = 1, 
                                                        aes(y= ..density.., fill = ..count..))+geom_density(fill="coral", alpha = 0.4)   + labs(x="", title = header)})
  
  ##boxplot--------------------------------------------------------------------------------------------------------
  output$boxplot <- renderPlot({
    dataset <- datasetInput()
    header <- datasetNames()
    Boxplot(y = dataset, main = header)
  })
  
  ##qqplot--------------------------------------------------------------------------------------------------------
  output$qqplot <- renderPlot({
    dataset <- datasetInput()
    header <- datasetNames()
    qqnorm(dataset[, 1], main = header) 
    qqline(dataset[, 1], col="coral")})
  
  ## correlation--------------------------------------------------------------------------------------------------------
  output$scatter <- renderPlot({
    pairs(swiss, lower.panel = panel.smooth, upper.panel = panel.cor,
          gap=0, row1attop=FALSE, main = "Scatterplot")})
  
  ##linear model-------------------------------------------------------------------------------------------------
  linearModel <- reactive({
    if (input$transformation == "Without Transformation") 
    {expln <- paste(input$checkbox, collapse = "+")
    as.formula(paste(input$regressand, "~", expln))} 
    else if (input$transformation == "Log(X)") 
    {expln <- paste("log(", input$checkbox, ")", collapse = "+")
    as.formula(paste(input$regressand, "~", expln))} 
    else if (input$transformation == "Log(Y)")
    {expln <- paste(input$checkbox, collapse = "+")
    as.formula(paste("log(",input$regressand, ")", "~", expln))}
    else if (input$transformation == "Log/Log")
    {expln <- paste("log(", input$checkbox, ")", collapse = "+")
    as.formula(paste("log(",input$regressand, ")", "~", expln))}
    else if (input$transformation == "Standardisation")
    {expln <- paste(input$checkbox, collapse = "+")
    as.formula(paste(input$regressand, "~", expln))} 
    else if (input$transformation == "X^2") {
      {expln <- paste(input$checkbox, "^2", collapse = "+")
      as.formula(paste(input$regressand, "~", expln))}
    }
    else if (input$transformation == "Y^2") {
      {expln <- paste(input$checkbox, collapse = "+")
      as.formula(paste(input$regressand, "^2", "~", expln))}
    }
    else if (input$transformation == "X^2/Y^2"){
      expln <- paste(input$checkbox, "^2", collapse = "+")
      as.formula(paste(input$regressand, "^2", "~", expln))}
  })
  
  ## @AIC ------------------------------
  output$AIC  <- renderPrint({
    linearModel()})
  
  output$stepmodel <- renderPrint({
    fit = lm(linearModel(), data=swiss[c(input$checkGroup),])
    step(fit)})
  
  ## @residue plots ------------------------------
  output$modelEquation <- renderPrint({
    linearModel()})
  
  output$residuePlot <- renderPlot ({
    if (input$transformation == "Standardisation"){
      fit = lm(linearModel(), data=normSwiss[c(input$checkGroup),]) 
      par(mfrow=c(2,2))
      plot(fit)
    }
    else{
      fit = lm(linearModel(), data=swiss[c(input$checkGroup),]) 
      par(mfrow=c(2,2))
      plot(fit)}
  })
  
  ## @model summary -----------------------------
  output$modelTheorem <- renderPrint({
    linearModel()})
  
  output$modelSummary <- renderPrint({
    if (input$transformation == "Standardisation"){
      fit = lm(linearModel(), data=normSwiss[c(input$checkGroup),]) 
      summary(fit)
    }
    else{
      fit = lm(linearModel(), data=swiss[c(input$checkGroup),]) 
      summary(fit)}
  })
  
  # ## @final model -----------------------------
  # output$finalModel <- renderText({
  #   fit = lm(linearModel(), data=swiss[c(input$checkGroup),])
  #   fit <- summary(fit)
  #   fit <- fit$coefficients
  #   alpha <- fit[1]
  #   beta <- fit[2]
  #   sprintf("y = %s + %s * x", alpha, beta)})
  
  output$modelPlot <- renderPrint({
    
    regressand <- as.data.frame(swiss[c(input$regressand)])
    cols <- length(swiss[, c(input$checkbox)])
    rows <- length(row.names(regressand))
    new_swiss <- matrix(ncol = cols, nrow = rows)
    rownames(new_swiss) <- row.names(regressand)
    colnames(new_swiss) <- colnames(swiss[c(input$checkbox),])
    for (i in 1:cols){
      tmp <- unlist(swiss[c(input$checkbox)][i], use.names = FALSE)
      new_swiss[, i] <- tmp
    }
    new_swiss[c(input$checkGroup),]
    # # plot(y = swiss[c(input$regressand)])
  })
  }

shinyApp(ui = ui, server = server)