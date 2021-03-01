library(shiny)
library(pheatmap)
library(shinyEventLogger)
library(Seurat)
library(ggplot2)
library(dplyr)
library(grid)
library(shinycssloaders)
library(shinythemes)
library(DT)
set_logging()
data <- readRDS("../treatments.so.rds")
data2 <- data; data3 <- data;
data4 <- data3; Idents(data4)<-data4@meta.data$Treatment
scaled <- readRDS("../my_data.rds")
# data2_scaled <- ScaleData(data2, vars.to.regress = c("batch_day","percent.mito","nUMI","percent.ribo"), 
#                           features = rownames(data2@assays$RNA@data))
Idents(data2)<-data2@meta.data$cellType_Treatment
Idents(data3)<-levels(data3@meta.data$Ind)
Idents(scaled)<-scaled@meta.data$cellType_Treatment
pairwise <- readRDS("../pairWiseTests.rds")

options(spinner.color="#74a5cf", spinner.color.background="#ffffff", spinner.size=1.5, spinner.type=4)

# Define UI for application that draws a histogram
ui <- navbarPage(theme = shinytheme("flatly"), title = tags$b("scRNA-seq Viz"), 
                 tabPanel("Documentation", value=-999,
                          mainPanel(imageOutput("docs")
                          )
                 ),
                 # Sidebar with a slider input for number of bins 
                 tabPanel(
                   "Data Clusters",
                   br(),
                   sidebarLayout(
                     sidebarPanel(
                       actionButton("refresh_1", "Refresh", class = "btn-success"),
                       selectInput("reduction", h3(tags$b("Reduction: ")), 
                                   choices = c("UMAP", "tSNE","PCA"), selected = "UMAP"),
                       selectInput("group.by", h3(tags$b("Group By: ")), 
                                   choices = c("cellType", "RNA_snn_res.2", "RNA_snn_res.1.5",
                                               "RNA_snn_res.1.25", "RNA_snn_res.1", "RNA_snn_res.0.8", "Donor", "Treatment", "Treatment vs. cellType"), selected = "cellType"),
                       br(),
                       sliderInput("pt.size", 
                                   label = "Point Size: ",
                                   min = 0.1, max = 1.5, value = 0.5),
                       br(),
                       downloadButton('downloadPlot1', 'Download Plot')
                     ),
                     
                     # Show a plot of the generated distribution
                     mainPanel(
                       withSpinner(plotOutput("DimPlot",  width = "100%"), hide.ui = FALSE)
                     )
                   )
                 ), 
                 tabPanel(
                   "Differentially Expressed Genes",
                   br(),
                   sidebarLayout(
                     sidebarPanel(
                       actionButton("refresh_2", "Refresh", class = "btn-success"),
                       conditionalPanel(
                         condition = "input.plot_type == 2",
                         h3(tags$b("Cell Type/Treatment Selection: ")),
                         br(),
                         h4(tags$b("Cell Types: ")),
                         br(),
                         checkboxInput("Basal_h", "Basal", value = TRUE),
                         checkboxInput("Ciliated_h", "Ciliated", value = FALSE),
                         checkboxInput("Intermediate_h", "Intermediate", value = FALSE),
                         checkboxInput("Preciliated_h", "Preciliated", value = FALSE),
                         checkboxInput("Proliferating_h", "Proliferating", value = FALSE),
                         checkboxInput("Secretory_h", "Secretory", value = FALSE),
                         br(),
                         selectInput("treatment_h", label = h4("Treatment: "), 
                                     choices = c("CO", "IFN", "IL13", "IL17"),
                                     selected = "CO"),
                         br(),
                         h3(tags$b("FDR/log2FC Cutoff: ")),
                         br(),
                         numericInput("FDR_h", label = "FDR cutoff: ", value = 0.1, min = 0, max = 1),
                         numericInput("log2FC_h", label = "log2FC cutoff: ", value = 1, min = -10, max = 10),
                         checkboxInput("include_downregulated_h", "Include only Down-regulated Genes", value = FALSE),
                         checkboxInput("include_upregulated_h", "Include only Up-regulated Genes", value = FALSE),
                         br(),
                         sliderInput("n_h", h3(tags$b("Number of genes: ")), min = 1, max = 100, value = 10),
                         br(),
                       ),
                       conditionalPanel(
                         condition = "input.plot_type != 2",
                         h3(tags$b("Cell Type/Treatment Selection: ")),
                         br(),
                         selectInput("cellType", label = "Cell Type: ", 
                                     choices = c("Basal", "Ciliated", "Intermediate", "Preciliated", "Proliferating_Basal", "Secretory"),
                                     selected = "Basal"),
                         selectInput("treatment", label = "Treatment: ", 
                                     choices = c("CO", "IFN", "IL13", "IL17"),
                                     selected = "CO"),
                         br(),
                         h3(tags$b("FDR/log2FC Cutoff: ")),
                         br(),
                         numericInput("FDR", label = "FDR cutoff: ", value = 0.1, min = 0, max = 1),
                         numericInput("log2FC", label = "log2FC cutoff: ", value = 1, min = -10, max = 10),
                         checkboxInput("include_downregulated", "Include only Down-regulated Genes", value = FALSE),
                         checkboxInput("include_upregulated", "Include only Up-regulated Genes", value = FALSE),
                         br(),
                         sliderInput("n", h3(tags$b("Number of genes: ")), min = 1, max = 100, value = 10),
                         br(),
                       ),
                       radioButtons("plot_type", h3(tags$b("Plot Type")),
                                    choices = list("Feature Plot" = 1, "Heat Map" = 2,
                                                   "Violin Plot" = 3, "Dot Plot" = 4), selected = 3),
                       br(),
                       sliderInput("pt.size2", 
                                   label = "Point Size: ",
                                   min = 0.1, max = 1.5, value = 0.5),
                       
                       
                       br(),
                       conditionalPanel(
                         condition = "input.plot_type != 2",
                         downloadButton('downloadPlot2', 'Download Plot')
                       ),
                       conditionalPanel(
                         condition = "input.plot_type == 2",
                         downloadButton('downloadPlot5', 'Download Plot'),
                         downloadButton('downloadData', 'Download Raw Data'),
                         br(),
                         checkboxInput("view_raw_h", "View Raw Data", value = FALSE)
                       )
                     ),
                     
                     # Show a plot of the generated distribution
                     mainPanel(
                       conditionalPanel(
                         condition = "input.plot_type == 1",
                         #plotOutput("FeaturePlot",  width = "100%"),
                         uiOutput("FeaturePlot.ui")
                       ),
                       conditionalPanel(
                         condition = "input.plot_type == 2",
                         conditionalPanel(
                           condition = "input.view_raw_h == 0",
                           withSpinner(plotOutput("HeatMap_h",  width = "100%"), hide.ui = FALSE)
                         ),
                         conditionalPanel(
                           condition = "input.view_raw_h == 1",
                           DT::dataTableOutput("table")
                         )
                       ),
                       conditionalPanel(
                         condition = "input.plot_type == 3",
                         #plotOutput("ViolinPlot",  width = "100%"),
                         uiOutput("ViolinPlot.ui")
                       ),
                       conditionalPanel(
                         condition = "input.plot_type == 4",
                         withSpinner(plotOutput("DotPlot",  width = "100%"), hide.ui = FALSE)
                       )
                     )
                   )
                 ),
                 tabPanel(
                   "Gene Analyzer",
                   br(),
                   sidebarLayout(
                     sidebarPanel(
                       actionButton("refresh_3", "Refresh", class = "btn-success"),
                       h3(tags$b("Gene Input: ")),
                       br(),
                       conditionalPanel(
                         condition = "input.plot_type2 == 4",
                         textInput("geneNameSpecial", label = "Gene Name: ", value = "POSTN")
                       ),
                       conditionalPanel(
                         condition = "input.plot_type2 == 5",
                         textInput("geneNameSpecial2", label = "Gene Name: ", value = "POSTN")
                       ),
                       conditionalPanel(
                         condition = "input.plot_type2 != 4 && input.plot_type2 != 5",
                         textInput("geneName", label = "Gene Name: ", value = "POSTN"),
                       ),
                       conditionalPanel(
                         condition = "input.plot_type2 == 5",
                         checkboxInput("make_log2FC", "Plot log2FC values", value = FALSE),
                       ),
                       br(),
                       h3(tags$b("Cell Type/Treatment Selection: ")),
                       br(),
                       conditionalPanel(
                         condition = "input.plot_type2 == 4 || input.plot_type2 == 5",
                         h3(tags$b("Cell Types: ")),
                         br(),
                         checkboxInput("Basal", "Basal", value = FALSE),
                         checkboxInput("Ciliated", "Ciliated", value = FALSE),
                         checkboxInput("Intermediate", "Intermediate", value = FALSE),
                         checkboxInput("Preciliated", "Preciliated", value = FALSE),
                         checkboxInput("Proliferating", "Proliferating", value = FALSE),
                         checkboxInput("Secretory", "Secretory", value = FALSE),
                         br(),
                         h3(tags$b("Treatments: ")),
                         br(),
                         checkboxInput("CO", "CO", value = FALSE),
                         checkboxInput("IFN", "IFN", value = FALSE),
                         checkboxInput("IL13", "IL13", value = FALSE),
                         checkboxInput("IL17", "IL17", value = FALSE),
                         checkboxInput("UT", "UT", value = TRUE),
                         br(),
                         conditionalPanel(
                           condition = "input.plot_type2 == 4",
                           h3(tags$b("Donor: ")),
                           br(),
                           checkboxInput("Sub1", "Donor 1", value = FALSE),
                           checkboxInput("Sub2", "Donor 2", value = FALSE),
                           checkboxInput("Sub3", "Donor 3", value = FALSE),
                           checkboxInput("Sub6", "Donor 4", value = FALSE),
                         ),
                       ),
                       conditionalPanel(
                         condition = "input.plot_type2 != 4 && input.plot_type2 != 5",
                         selectInput("cellType2", label = "Cell Type: ", 
                                     choices = c("Basal", "Ciliated", "Intermediate", "Preciliated", "Proliferating_Basal", "Secretory"),
                                     selected = "Basal"),
                         selectInput("treatment2", label = "Treatment: ", 
                                     choices = c("CO", "IFN", "IL13", "IL17"),
                                     selected = "CO"),
                       ),
                       br(),
                       radioButtons("plot_type2", h3(tags$b("Plot Type")),
                                    choices = list("Feature Plot" = 1, "Violin Plot (Cell Type Expression)" = 2,
                                                   "Violin Plot (Untreated/Treated)" = 3, "Dot Plot" = 4, "HeatMap" = 5), selected = 1),
                       br(),
                       sliderInput("pt.size3", 
                                   label = "Point Size: ",
                                   min = 0.1, max = 1.5, value = 0.5),
                       br(),
                       conditionalPanel(
                         condition = "input.plot_type2 == 4",
                         downloadButton('downloadPlot4', 'Download Plot'),
                         downloadButton('downloadData3', 'Download Raw Data'),
                         br(),
                         checkboxInput("view_raw_d", "View Raw Data", value = FALSE)
                       ),
                       conditionalPanel(
                         condition = "input.plot_type2 != 4",
                         downloadButton('downloadPlot3', 'Download Plot')
                       )
                     ),
                     
                     # Show a plot of the generated distribution
                     mainPanel(
                       conditionalPanel(
                         condition = "input.plot_type2 == 1",
                         #plotOutput("FeaturePlot2",  width = "100%"),
                         uiOutput("FeaturePlot2.ui")
                       ),
                       conditionalPanel(
                         condition = "input.plot_type2 == 2",
                         #plotOutput("ViolinPlot2CellTypeExpression",  width = "100%"),
                         uiOutput("ViolinPlot2CellTypeExpression.ui")
                       ),
                       conditionalPanel(
                         condition = "input.plot_type2 == 3",
                         #plotOutput("ViolinPlot2UntreatedTreated",  width = "100%"),
                         uiOutput("ViolinPlot2UntreatedTreated.ui")
                       ),
                       conditionalPanel(
                         condition = "input.plot_type2 == 4",
                         conditionalPanel(
                           condition = "input.view_raw_d == 0",
                           withSpinner(plotOutput("DotPlotSpecial",  width = "100%"), hide.ui = FALSE)
                         ),
                         conditionalPanel(
                           condition = "input.view_raw_d == 1",
                           DT::dataTableOutput("table2")
                         )
                       ),
                       conditionalPanel(
                         condition = "input.plot_type2 == 5",
                         withSpinner(plotOutput("HeatMapSpecial",  width = "100%"), hide.ui = FALSE)
                       ),
                       br(),
                       textOutput("validation"),
                       br()
                     )
                   )
                 ),
                 tabPanel(
                   "Marker Genes",
                   br(),
                   sidebarLayout(
                     sidebarPanel(
                       actionButton("refresh_4", "Refresh", class = "btn-success"),
                       h3(tags$b("Cell Type/Treatment Selection: ")),
                       br(),
                       h4(tags$b("Cell Types: ")),
                       br(),
                       checkboxInput("Basal_m", "Basal", value = TRUE),
                       checkboxInput("Ciliated_m", "Ciliated", value = FALSE),
                       checkboxInput("Intermediate_m", "Intermediate", value = FALSE),
                       checkboxInput("Preciliated_m", "Preciliated", value = FALSE),
                       checkboxInput("Proliferating_m", "Proliferating", value = FALSE),
                       checkboxInput("Secretory_m", "Secretory", value = FALSE),
                       br(),
                       selectInput("treatment_m", label = h4("Treatment: "), 
                                   choices = c("CO", "IFN", "IL13", "IL17", "UT"),
                                   selected = "CO"),
                       br(),
                       checkboxInput("include_downregulated_m", "Include only Down-regulated Genes", value = FALSE),
                       checkboxInput("include_upregulated_m", "Include only Up-regulated Genes", value = FALSE),
                       br(),
                       sliderInput("n_m", h3(tags$b("Number of genes: ")), min = 1, max = 100, value = 10),
                       br(),
                       downloadButton('downloadPlot6', 'Download Plot'),
                       downloadButton('downloadData2', 'Download Raw Data'),
                       br(),
                       checkboxInput("view_raw_m", "View Raw Data", value = FALSE)
                     ),
                     
                     # Show a plot of the generated distribution
                     mainPanel(
                       conditionalPanel(
                         condition = "input.view_raw_m == 0",
                         withSpinner(plotOutput("DotPlot_m",  width = "100%"), hide.ui = FALSE)
                       ),
                       conditionalPanel(
                         condition = "input.view_raw_m == 1",
                         DT::dataTableOutput("table3")
                       )
                     )
                   )
                 )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  set_logging_session() 
  log_event("Hello World!")
  plotCount <- reactive({
    input$refresh_2
    as.numeric(isolate(input$n))
  })
  value <- reactiveVal(0) 
  
  plotHeight <- reactive(250*ceiling((plotCount())/2))   
  plotHeight2 <- reactive(500*ceiling((value())/2))   
  
  output$DotPlotSpecial <- renderPlot({
    input$refresh_3
    # data is what reduction you are using
    gene <- strsplit(isolate(input$geneNameSpecial), ", ")
    SUB1<-subset(data3, idents=c("Sub1")); SUB2<-subset(data3, idents=c("Sub2")); 
    SUB3<-subset(data3, idents=c("Sub3")); SUB6<-subset(data3, idents=c("Sub6"));
    Sub1 <- isolate(input$Sub1); Sub2 <- isolate(input$Sub2); Sub3 <- isolate(input$Sub3); Sub6 <- isolate(input$Sub6);
    
    types <- c(); treats <- c()
    if (isolate(input$Basal)){types <- append(types, "Basal")}
    if (isolate(input$Ciliated)){types <- append(types, "Ciliated")}
    if (isolate(input$Intermediate)){types <- append(types, "Intermediate")}
    if (isolate(input$Preciliated)){types <- append(types, "Preciliated")}
    if (isolate(input$Proliferating)){types <- append(types, "Proliferating")}
    if (isolate(input$Secretory)){types <- append(types, "Secretory")}
    if (isolate(input$CO)){treats <- append(treats, "CO")}
    if (isolate(input$IFN)){treats <- append(treats, "IFN")}
    if (isolate(input$IL13)){treats <- append(treats, "IL13")}
    if (isolate(input$IL17)){treats <- append(treats, "IL17")}
    if (isolate(input$UT)){treats <- append(treats, "UT")}
    if (length(types) == 0){
      if (length(treats) == 0){
        ret = 0; Idents(SUB1) <- levels(data); Idents(SUB2) <- levels(data); Idents(SUB3) <- levels(data); Idents(SUB6) <- levels(data);   
        if (Sub1 && length(ret)== 1){ret <- DotPlot(SUB1, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub1")}
        else if (Sub1){ret = ret + DotPlot(SUB1, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub1")}
        if (Sub2 && length(ret)== 1){ret <- DotPlot(SUB2, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub2")}
        else if (Sub2){ret = ret + DotPlot(SUB2, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub2")}
        if (Sub3 && length(ret)== 1){ret <- DotPlot(SUB3, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub3")}
        else if (Sub3){ret = ret + DotPlot(SUB3, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub3")}
        if (Sub6 && length(ret)== 1){ret <- DotPlot(SUB6, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub6")}
        else if (Sub6){ret = ret + DotPlot(SUB6, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub6")}
        if (!Sub6 && !Sub2 && !Sub3 && !Sub1){DotPlot(data, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))}
        else {ret}
      }else {
        ret = 0; Idents(SUB1) <- levels(data4); Idents(SUB2) <- levels(data4); Idents(SUB3) <- levels(data4); Idents(SUB6) <- levels(data4); 
        if (Sub1 && length(ret)== 1){ret <- DotPlot(SUB1, features = gene[[1]], idents = treats, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub1")}
        else if (Sub1){ret = ret + DotPlot(SUB1, features = gene[[1]], idents = treats, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub1")}
        if (Sub2 && length(ret)== 1){ret <- DotPlot(SUB2, features = gene[[1]], idents = treats, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub2")}
        else if (Sub2){ret = ret + DotPlot(SUB2, features = gene[[1]], idents = treats, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub2")}
        if (Sub3 && length(ret)== 1){ret <- DotPlot(SUB3, features = gene[[1]], idents = treats, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub3")}
        else if (Sub3){ret = ret + DotPlot(SUB3, features = gene[[1]], idents = treats, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub3")}
        if (Sub6 && length(ret)== 1){ret <- DotPlot(SUB6, features = gene[[1]], idents = treats, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub6")}
        else if (Sub6){ret = ret + DotPlot(SUB6, features = gene[[1]], idents = treats, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub6")}
        if (!Sub6 && !Sub2 && !Sub3 && !Sub1){DotPlot(data4, features = gene[[1]], idents =treats, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))}
        else {ret}
      }
    }else if (length(treats) == 0){
      if (length(types) == 0){
        ret = 0; Idents(SUB1) <- levels(data); Idents(SUB2) <- levels(data); Idents(SUB3) <- levels(data); Idents(SUB6) <- levels(data);   
        if (Sub1 && length(ret)== 1){ret <- DotPlot(SUB1, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub1")}
        else if (Sub1){ret = ret + DotPlot(SUB1, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub1")}
        if (Sub2 && length(ret)== 1){ret <- DotPlot(SUB2, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub2")}
        else if (Sub2){ret = ret + DotPlot(SUB2, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub2")}
        if (Sub3 && length(ret)== 1){ret <- DotPlot(SUB3, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub3")}
        else if (Sub3){ret = ret + DotPlot(SUB3, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub3")}
        if (Sub6 && length(ret)== 1){ret <- DotPlot(SUB6, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub6")}
        else if (Sub6){ret = ret + DotPlot(SUB6, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub6")}
        if (!Sub6 && !Sub2 && !Sub3 && !Sub1){DotPlot(data, features = gene[[1]], idents = levels(data), cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))}
        else {ret}
      }else {
        ret = 0; Idents(SUB1) <- levels(data); Idents(SUB2) <- levels(data); Idents(SUB3) <- levels(data); Idents(SUB6) <- levels(data);   
        if (Sub1 && length(ret)== 1){ret <- DotPlot(SUB1, features = gene[[1]], idents = types, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub1")}
        else if (Sub1){ret = ret + DotPlot(SUB1, features = gene[[1]], idents = types, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub1")}
        if (Sub2 && length(ret)== 1){ret <- DotPlot(SUB2, features = gene[[1]], idents = types, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub2")}
        else if (Sub2){ret = ret + DotPlot(SUB2, features = gene[[1]], idents = types, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub2")}
        if (Sub3 && length(ret)== 1){ret <- DotPlot(SUB3, features = gene[[1]], idents = types, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub3")}
        else if (Sub3){ret = ret + DotPlot(SUB3, features = gene[[1]], idents = types, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub3")}
        if (Sub6 && length(ret)== 1){ret <- DotPlot(SUB6, features = gene[[1]], idents = types, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub6")}
        else if (Sub6){ret = ret + DotPlot(SUB6, features = gene[[1]], idents = types, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub6")}
        if (!Sub6 && !Sub2 && !Sub3 && !Sub1){DotPlot(data, features = gene[[1]], idents = types, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))}
        else {ret}
      }
    }else {
      Idents(data2) <- sort(levels(Idents(data2)))
      ret = 0; Idents(SUB1) <- levels(data2); Idents(SUB2) <- levels(data2); Idents(SUB3) <- levels(data2); Idents(SUB6) <- levels(data2);
      combo <- c()
      for (type in types){
        for (treat in treats){
          combo <- append(combo, paste(type, treat, sep="_"))
        }
      }
      if (Sub1 && length(ret)== 1){ret <- DotPlot(SUB1, features = gene[[1]], idents = combo, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub1")}
      else if (Sub1){ret = ret + DotPlot(SUB1, features = gene[[1]], idents = combo, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub1")}
      if (Sub2 && length(ret)== 1){ret <- DotPlot(SUB2, features = gene[[1]], idents = combo, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub2")}
      else if (Sub2){ret = ret + DotPlot(SUB2, features = gene[[1]], idents = combo, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub2")}
      if (Sub3 && length(ret)== 1){ret <- DotPlot(SUB3, features = gene[[1]], idents = combo, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub3")}
      else if (Sub3){ret = ret + DotPlot(SUB3, features = gene[[1]], idents = combo, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub3")}
      if (Sub6 && length(ret)== 1){ret <- DotPlot(SUB6, features = gene[[1]], idents = combo, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub6")}
      else if (Sub6){ret = ret + DotPlot(SUB6, features = gene[[1]], idents = combo, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Sub6")}
      if (!Sub6 && !Sub2 && !Sub3 && !Sub1){DotPlot(data2, features = gene[[1]], idents = combo, cols = c("red", "blue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))}
      else {ret}
    }
    
  }, height = 1000, width = 900)
  
  output$HeatMapSpecial <- renderPlot({
    input$refresh_3
    # data is what reduction you are using
    gene <- strsplit(isolate(input$geneNameSpecial2), ", ")
    types <- c(); treats <- c()
    flag <- isolate(input$make_log2FC)
    if (isolate(input$Basal)){types <- append(types, "Basal")}
    if (isolate(input$Ciliated)){types <- append(types, "Ciliated")}
    if (isolate(input$Intermediate)){types <- append(types, "Intermediate")}
    if (isolate(input$Preciliated)){types <- append(types, "Preciliated")}
    if (isolate(input$Proliferating)){types <- append(types, "Proliferating")}
    if (isolate(input$Secretory)){types <- append(types, "Secretory")}
    if (isolate(input$CO)){treats <- append(treats, "CO")}
    if (isolate(input$IFN)){treats <- append(treats, "IFN")}
    if (isolate(input$IL13)){treats <- append(treats, "IL13")}
    if (isolate(input$IL17)){treats <- append(treats, "IL17")}
    if (isolate(input$UT)){treats <- append(treats, "UT")}
    validate(
      need((length(types) != 0 || length(treats) != 0), "Must select a treatment or cell-type")
    )
    validate(
      need(!(flag && (length(treats) != 1 || length(types) == 0)), "Must select a one treatment AND cell-type combination for log2FC analysis")
    )
    log_event(flag)
    if (flag){
      log_event("log2FC: YES!")
      run <- 0;
      mat <- matrix(1:(length(gene[[1]]) * length(types)), nrow = length(gene[[1]]), ncol = length(types))
      for (g in c(1: length(gene[[1]]))){
        for (type in types){
          treatment_name <- paste(paste(type, treats[1], sep="_"), paste(type, "UT", sep="_"), sep=".vs.")
          val <- pairwise[[treatment_name]]
          ind <- -1
          for (i in c(1: nrow(val))){
            if (val$Gene[i] == gene[[1]][g]){
              ind <- i
              break;
            }
          }
          if (ind == -1){
            mat[g, run + 1] <- 0
          }else {
            mat[g, run + 1] <- val$log2FC[ind]
          }
          run <- (run + 1) %% length(types)
        }
      }
      rownames(mat) = gene[[1]]
      colnames(mat) = types
      breaksList = seq(-1, 2, by = 0.01)
      pheatmap(mat, breaks = breaksList, color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)), cluster_rows = FALSE, cluster_cols = FALSE)
    }else {
      if (length(types) == 0){
        Idents(scaled)<-scaled@meta.data$Treatment
        ThreeCellTypes <- subset(scaled, idents = treats)
        levels(ThreeCellTypes) <- sort(levels(ThreeCellTypes))
        cluster.averages <- AverageExpression(ThreeCellTypes, return.seurat = TRUE)
        DoHeatmap(cluster.averages, features = gene[[1]], draw.lines = FALSE) + scale_fill_gradientn(colors = c("blue", "white", "red"))
      }else if (length(treats) == 0){
        Idents(scaled)<-scaled@meta.data$cellType
        ThreeCellTypes <- subset(scaled, idents = types)
        levels(ThreeCellTypes) <- sort(levels(ThreeCellTypes))
        cluster.averages <- AverageExpression(ThreeCellTypes, return.seurat = TRUE)
        DoHeatmap(cluster.averages, features = gene[[1]], draw.lines = FALSE) + scale_fill_gradientn(colors = c("blue", "white", "red"))
      }else {
        combo <- c()
        for (type in types){
          for (treat in treats){
            combo <- append(combo, paste(type, treat, sep="_"))
          }
        }
        Idents(scaled)<-scaled@meta.data$cellType_Treatment
        ThreeCellTypes <- subset(scaled, idents = combo)
        levels(ThreeCellTypes) <- sort(levels(ThreeCellTypes))
        cluster.averages <- AverageExpression(ThreeCellTypes, return.seurat = TRUE)
        DoHeatmap(cluster.averages, features = gene[[1]], draw.lines = FALSE) + scale_fill_gradientn(colors = c("blue", "white", "red"))
      }
    }
  }, height = 1000, width = 900)
  
  output$DimPlot <- renderPlot({
    # data is what reduction you are using
    input$refresh_1
    reduc <- switch(isolate(input$reduction), 
                    "UMAP" = "umap",
                    "tSNE" = "tsne",
                    "PCA" = "pca")
    group.by <- isolate(input$group.by)
    pt.size <- isolate(input$pt.size)
    if (group.by == "Donor"){
      group.by = "Ind"
      DimPlot(object = data, group.by=group.by, pt.size=pt.size, reduction = reduc, label = T)
    }else if (group.by == "Treatment vs. cellType"){
      group.by = "cellType_Treatment"
      ThreeCellTypes<-subset(data, idents=c("Basal", "Ciliated", "Secretory"))
      Idents(ThreeCellTypes) <- ThreeCellTypes@meta.data$cellType_Treatment
      DimPlot(object = ThreeCellTypes, group.by=group.by, pt.size=pt.size, reduction = reduc, label = T)
    }else {
      DimPlot(object = data, group.by=group.by, pt.size=pt.size, reduction = reduc, label = T)
    }
  }, height = 700, width = 900)
  
  output$DotPlot <- renderPlot({
    input$refresh_2
    # data is what reduction you are using
    cellType <- isolate(input$cellType)
    treatment <- isolate(input$treatment)
    fdr <- isolate(input$FDR)
    n <- isolate(input$n)
    log2FC <- isolate(input$log2FC)
    treatment <- isolate(input$treatment)
    pt.size <- isolate(input$pt.size2)
    include_downregulated <- isolate(input$include_downregulated)
    include_upregulated <- isolate(input$include_upregulated)
    treatment_name <- paste(paste(cellType, treatment, sep="_"),
                            paste(cellType, "UT", sep="_"), sep=".vs.")
    getSigGenes_upregulated <- function(FDRcutoff, log2fcutoff, dataframeOfCoef = df){
      validate(
        need((include_downregulated != TRUE || include_upregulated != TRUE), "Cannot plot only downregulated and only upregulated genes at the same time")
      )
      if (include_downregulated == TRUE){
        newdata <- df[order(df$log2FC), ]
        dfsub <- subset(newdata, ((newdata$log2FC < 0)&(abs(newdata$log2FC) > log2fcutoff)&(newdata$hurdle_fdr<FDRcutoff)))
      }else if (include_upregulated == TRUE){
        newdata <- df[order(-df$log2FC), ]
        dfsub <- subset(newdata, ((newdata$log2FC > log2fcutoff)&(newdata$hurdle_fdr<FDRcutoff)))
      }else {
        newdata <- df[order(-abs(df$log2FC)), ]
        dfsub <- subset(newdata, ((abs(newdata$log2FC) > log2fcutoff)&(newdata$hurdle_fdr<FDRcutoff)))
      }
      return (dfsub[, 1])
    }
    df <- pairwise[[treatment_name]]
    top_n <- getSigGenes_upregulated(fdr, log2FC, df)
    top_n <- top_n[1: n]
    validate(
      need(top_n[1] != "NA", "It looks like no differentially expressed genes were found, consider changing the FDR/log2FC cutoffs")
    )
    # draw the FeaturePlot
    DotPlot(data2, features = top_n, idents = c(paste(cellType, treatment, sep="_"), paste(cellType, "UT", sep="_"))) + coord_flip()
  }, height = 1000, width = 900)
  
  output$FeaturePlot <- renderPlot({
    input$refresh_2
    # data is what reduction you are using
    cellType <- isolate(input$cellType)
    treatment <- isolate(input$treatment)
    fdr <- isolate(input$FDR)
    n <- isolate(input$n)
    log2FC <- isolate(input$log2FC)
    treatment <- isolate(input$treatment)
    pt.size <- isolate(input$pt.size2)
    include_downregulated <- isolate(input$include_downregulated)
    include_upregulated <- isolate(input$include_upregulated)
    treatment_name <- paste(paste(cellType, treatment, sep="_"),
                            paste(cellType, "UT", sep="_"), sep=".vs.")
    getSigGenes_upregulated <- function(FDRcutoff, log2fcutoff, dataframeOfCoef = df){
      validate(
        need((include_downregulated != TRUE || include_upregulated != TRUE), "Cannot plot only downregulated and only upregulated genes at the same time")
      )
      if (include_downregulated == TRUE){
        newdata <- df[order(df$log2FC), ]
        dfsub <- subset(newdata, ((newdata$log2FC < 0)&(abs(newdata$log2FC) > log2fcutoff)&(newdata$hurdle_fdr<FDRcutoff)))
      }else if (include_upregulated == TRUE){
        newdata <- df[order(-df$log2FC), ]
        dfsub <- subset(newdata, ((newdata$log2FC > log2fcutoff)&(newdata$hurdle_fdr<FDRcutoff)))
      }else {
        newdata <- df[order(-abs(df$log2FC)), ]
        dfsub <- subset(newdata, ((abs(newdata$log2FC) > log2fcutoff)&(newdata$hurdle_fdr<FDRcutoff)))
      }
      return (dfsub[, 1])
    }
    df <- pairwise[[treatment_name]]
    top_n <- getSigGenes_upregulated(fdr, log2FC, df)
    top_n <- top_n[1: n]
    validate(
      need(top_n[1] != "NA", "It looks like no differentially expressed genes were found, consider changing the FDR/log2FC cutoffs")
    )
    # draw the FeaturePlot
    FeaturePlot(data, features = top_n, pt.size=pt.size)
  }, width = 900)
  output$FeaturePlot.ui <- renderUI({
    withSpinner(plotOutput("FeaturePlot", height = plotHeight()), hide.ui = FALSE)
  })
  
  output$FeaturePlot2 <- renderPlot({
    input$refresh_3
    # data is what reduction you are using
    gene <- strsplit(isolate(input$geneName), ", ")
    pt.size <- isolate(input$pt.size3)
    # draw the FeaturePlot
    value(length(gene[[1]]))
    FeaturePlot(data, features = gene[[1]], pt.size=pt.size)
  }, width = 900)
  output$FeaturePlot2.ui <- renderUI({
    input$refresh_3
    gene <- strsplit(isolate(input$geneName), ", ")
    validate(
      need(length(gene[[1]]) != "0", "Please input a gene name")
    )
    p <- sort(data@assays$RNA@var.features)
    num <- 0
    for (i in gene[[1]]){l <- 0; r <- 10000
    while (l < r){mi <- round((l + r)/2, 0); if (p[mi] < i){l = mi + 1}else {r = mi - 1}}
    if (p[l] != i && p[l - 1] != i && p[l + 1] != i){num <- num + 1}}
    value(length(gene[[1]]) - num)
    withSpinner(plotOutput("FeaturePlot2", height = plotHeight2()), hide.ui = FALSE)
  })
  
  output$validation <- renderText({
    input$refresh_3
    # data is what reduction you are using
    gene <- strsplit(isolate(input$geneName), ", ")
    ret <- ""; p <- sort(data@assays$RNA@var.features)
    for (i in gene[[1]]){l <- 0; r <- 10000
      while (l < r){mi <- round((l + r)/2, 0); if (p[mi] < i){l = mi + 1}else {r = mi - 1}}
      if (p[l] != i && p[l - 1] != i && p[l + 1] != i){if (nchar(ret) != 0){ret <- paste(ret, i, sep = ", ")}else{ret <- i}}}
    if (nchar(ret) != 0){ret <- paste("These genes were not found:", ret); }
    ret
  })
  
  output$ViolinPlot2UntreatedTreated <- renderPlot({
    input$refresh_3
    # data is what reduction you are using
    gene <- strsplit(isolate(input$geneName), ", ")
    cellType <- isolate(input$cellType2)
    treatment <- isolate(input$treatment2)
    pt.size <- isolate(input$pt.size3)
    treatment_name <- paste(paste(cellType, treatment, sep="_"),
                            paste(cellType, "UT", sep="_"), sep=".vs.")
    # draw the FeaturePlot
    value(length(gene[[1]]))
    VlnPlot(data2, features = gene[[1]], pt.size=pt.size, idents = c(paste(cellType, treatment, sep="_"), paste(cellType, "UT", sep="_")))
  }, width = 900)
  output$ViolinPlot2UntreatedTreated.ui <- renderUI({
    input$refresh_3
    gene <- strsplit(isolate(input$geneName), ", ")
    validate(
      need(length(gene[[1]]) != "0", "Please input a gene name")
    )
    p <- sort(data@assays$RNA@var.features)
    num <- 0
    for (i in gene[[1]]){l <- 0; r <- 10000
    while (l < r){mi <- round((l + r)/2, 0); if (p[mi] < i){l = mi + 1}else {r = mi - 1}}
    if (p[l] != i && p[l - 1] != i && p[l + 1] != i){num <- num + 1}}
    value(length(gene[[1]]) - num)
    withSpinner(plotOutput("ViolinPlot2UntreatedTreated", height = plotHeight2()), hide.ui = FALSE)
  })
  
  output$ViolinPlot2CellTypeExpression <- renderPlot({
    input$refresh_3
    # data is what reduction you are using
    gene <- strsplit(isolate(input$geneName), ", ")
    pt.size <- isolate(input$pt.size3)
    value(length(gene[[1]]))
    VlnPlot(data, features = gene[[1]], pt.size=pt.size)
  }, width = 900)
  output$ViolinPlot2CellTypeExpression.ui <- renderUI({
    input$refresh_3
    gene <- strsplit(isolate(input$geneName), ", ")
    validate(
      need(length(gene[[1]]) != "0", "Please input a gene name")
    )
    p <- sort(data@assays$RNA@var.features)
    num <- 0
    for (i in gene[[1]]){l <- 0; r <- 10000
    while (l < r){mi <- round((l + r)/2, 0); if (p[mi] < i){l = mi + 1}else {r = mi - 1}}
    if (p[l] != i && p[l - 1] != i && p[l + 1] != i){num <- num + 1}}
    value(length(gene[[1]]) - num)
    withSpinner(plotOutput("ViolinPlot2CellTypeExpression", height = plotHeight2()), hide.ui = FALSE)
  })
  
  output$ViolinPlot <- renderPlot({
    input$refresh_2
    # data is what reduction you are using
    cellType <- isolate(input$cellType)
    treatment <- isolate(input$treatment)
    fdr <- isolate(input$FDR)
    n <- isolate(input$n)
    log2FC <- isolate(input$log2FC)
    treatment <- isolate(input$treatment)
    pt.size <- isolate(input$pt.size2)
    include_downregulated <- isolate(input$include_downregulated)
    include_upregulated <- isolate(input$include_upregulated)
    treatment_name <- paste(paste(cellType, treatment, sep="_"),
                            paste(cellType, "UT", sep="_"), sep=".vs.")
    getSigGenes_upregulated <- function(FDRcutoff, log2fcutoff, dataframeOfCoef = df){
      validate(
        need((include_downregulated != TRUE || include_upregulated != TRUE), "Cannot plot only downregulated and only upregulated genes at the same time")
      )
      if (include_downregulated == TRUE){
        newdata <- df[order(df$log2FC), ]
        dfsub <- subset(newdata, ((newdata$log2FC < 0)&(abs(newdata$log2FC) > log2fcutoff)&(newdata$hurdle_fdr<FDRcutoff)))
      }else if (include_upregulated == TRUE){
        newdata <- df[order(-df$log2FC), ]
        dfsub <- subset(newdata, ((newdata$log2FC > log2fcutoff)&(newdata$hurdle_fdr<FDRcutoff)))
      }else {
        newdata <- df[order(-abs(df$log2FC)), ]
        dfsub <- subset(newdata, ((abs(newdata$log2FC) > log2fcutoff)&(newdata$hurdle_fdr<FDRcutoff)))
      }
      return (dfsub[, 1])
    }
    df <- pairwise[[treatment_name]]
    top_n <- getSigGenes_upregulated(fdr, log2FC, df)
    top_n <- top_n[1: n]
    validate(
      need(top_n[1] != "NA", "It looks like no differentially expressed genes were found, consider changing the FDR/log2FC cutoffs")
    )
    # draw the FeaturePlot
    p <- VlnPlot(data2, features = top_n, pt.size=pt.size, idents = c(paste(cellType, treatment, sep="_"), paste(cellType, "UT", sep="_")), combine = FALSE, sort = TRUE)
    avgExpression <- function(geneName, Treatment, data){
      avg <- AverageExpression(object=data, features = geneName, slot = "scale.data", verbose= "FALSE")
      return (avg$RNA[[Treatment]])
    }
    media <- function(data){
      if (length(data) %% 2 == 0){
        return ((data[length(data)/2] + data[length(data)/2 + 1])/2)
      }
      return (data[ceiling(length(data)/2)])
    }
    for (i in 1:length(p)){
      #Idents(data2)<-c(paste(cellType, treatment, sep="_"))
      untreated <- c();
      treated <- c();
      for (j in c(1:length(p[[i]]$data[[1]]))){
        c <- levels(droplevels(p[[i]]$data[[2]][j]))
        len <- nchar(c); s <- substring(c,len-1,len)
        if (s == "UT"){
          untreated <- append(untreated, p[[i]]$data[[1]][j]);
        }else {
          treated <- append(treated, p[[i]]$data[[1]][j]);
        }
      }
      sort(untreated); sort(treated)
      placeholder <- signif(media(treated),digits=4)
      placeholder2 <- signif(media(untreated),digits=4)
      if (placeholder < placeholder2){
        tmp <- placeholder
        placeholder <- placeholder2
        placeholder2 <- tmp
      }
      p[[i]] = p[[i]] + labs(fill = "Median Expression") + scale_fill_discrete(name = "Median Expression", labels = c(placeholder, placeholder2))
    }
    #Idents(data2)<-data2@meta.data$cellType_Treatment
    if (length(p) > 1){
      for (i in 1:length(p)){
        if (i %% 2 == 0){
          p[[i - 1]] = p[[i - 1]] | p[[i]]
        }
      }
      if (length(p) %% 2 != 0){
        p[[length(p)]] = p[[length(p)]] | (ggplot() + theme_void())
      }
      for (i in 2:length(p)){
        if (i %% 2 == 1){
          p[[1]] = p[[1]] / p[[i]]
        }
      }
    }
    p[[1]]
  }, width = 900)
  output$ViolinPlot.ui <- renderUI({
    withSpinner(plotOutput("ViolinPlot", height = plotHeight()), hide.ui = FALSE)
  })
  
  
  output$HeatMap <- renderPlot({
    input$refresh_2
    # data is what reduction you are using
    cellType <- isolate(input$cellType)
    treatment <- isolate(input$treatment)
    fdr <- isolate(input$FDR)
    n <- isolate(input$n)
    log2FC <- isolate(input$log2FC)
    treatment <- isolate(input$treatment)
    include_downregulated <- isolate(input$include_downregulated)
    include_upregulated <- isolate(input$include_upregulated)
    treatment_name <- paste(paste(cellType, treatment, sep="_"),
                            paste(cellType, "UT", sep="_"), sep=".vs.")
    getSigGenes_upregulated <- function(FDRcutoff, log2fcutoff, dataframeOfCoef = df){
      validate(
        need((include_downregulated != TRUE || include_upregulated != TRUE), "Cannot plot only downregulated and only upregulated genes at the same time")
      )
      if (include_downregulated == TRUE){
        newdata <- df[order(df$log2FC), ]
        dfsub <- subset(newdata, ((newdata$log2FC < 0)&(abs(newdata$log2FC) > log2fcutoff)&(newdata$hurdle_fdr<FDRcutoff)))
      }else if (include_upregulated == TRUE){
        newdata <- df[order(-df$log2FC), ]
        dfsub <- subset(newdata, ((newdata$log2FC > log2fcutoff)&(newdata$hurdle_fdr<FDRcutoff)))
      }else {
        newdata <- df[order(-abs(df$log2FC)), ]
        dfsub <- subset(newdata, ((abs(newdata$log2FC) > log2fcutoff)&(newdata$hurdle_fdr<FDRcutoff)))
      }
      return (dfsub[, 1])
    }
    df <- pairwise[[treatment_name]]
    top_n <- getSigGenes_upregulated(fdr, log2FC, df)
    top_n <- top_n[1: n]
    validate(
      need(top_n[1] != "NA", "It looks like no differentially expressed genes were found, consider changing the FDR/log2FC cutoffs")
    )
    # draw the FeaturePlot
    ThreeCellTypes <- subset(scaled, idents = c(paste(cellType, treatment, sep="_"), paste(cellType, "UT", sep="_")))
    DoHeatmap(ThreeCellTypes, features=top_n, draw.lines = FALSE)
  }, height = 1000, width = 900)
  
  output$HeatMap_h <- renderPlot({
    input$refresh_2
    # data is what reduction you are using
    types <- c()
    if (isolate(input$Basal_h)){types <- append(types, "Basal")}
    if (isolate(input$Ciliated_h)){types <- append(types, "Ciliated")}
    if (isolate(input$Intermediate_h)){types <- append(types, "Intermediate")}
    if (isolate(input$Preciliated_h)){types <- append(types, "Preciliated")}
    if (isolate(input$Proliferating_h)){types <- append(types, "Proliferating")}
    if (isolate(input$Secretory_h)){types <- append(types, "Secretory")}
    validate(
      need(length(types) != 0, "Must have at least one cell type")
    )
    treatment <- isolate(input$treatment_h)
    fdr <- isolate(input$FDR_h)
    n <- isolate(input$n_h)
    log2FC <- isolate(input$log2FC_h)
    include_downregulated <- isolate(input$include_downregulated_h)
    include_upregulated <- isolate(input$include_upregulated_h)
    getSigGenes_upregulated <- function(FDRcutoff, log2fcutoff, dataframeOfCoef = df){
      validate(
        need((include_downregulated != TRUE || include_upregulated != TRUE), "Cannot plot only downregulated and only upregulated genes at the same time")
      )
      if (include_downregulated == TRUE){
        newdata <- df[order(df$log2FC), ]
        dfsub <- subset(newdata, ((newdata$log2FC < 0)&(abs(newdata$log2FC) > log2fcutoff)&(newdata$hurdle_fdr<FDRcutoff)))
      }else if (include_upregulated == TRUE){
        newdata <- df[order(-df$log2FC), ]
        dfsub <- subset(newdata, ((newdata$log2FC > log2fcutoff)&(newdata$hurdle_fdr<FDRcutoff)))
      }else {
        newdata <- df[order(-abs(df$log2FC)), ]
        dfsub <- subset(newdata, ((abs(newdata$log2FC) > log2fcutoff)&(newdata$hurdle_fdr<FDRcutoff)))
      }
      return (dfsub[, 1])
    }
    top_n <- c(); idents <- c();
    for (i in 1:length(types)){
      treatment_name <- paste(paste(types[i], treatment, sep="_"), paste(types[i], "UT", sep="_"), sep=".vs.")
      df <- pairwise[[treatment_name]]
      top_n_curr <- getSigGenes_upregulated(fdr, log2FC, df)
      top_n_curr <- top_n_curr[1: n]
      top_n <- union(top_n, top_n_curr)
      idents <- append(idents, paste(types[i], treatment, sep="_")); idents <- append(idents, paste(types[i], "UT", sep="_"))
    }
    validate(
      need(top_n[1] != "NA", "It looks like no differentially expressed genes were found, consider changing the FDR/log2FC cutoffs")
    )
    # draw the FeaturePlot
    ThreeCellTypes <- subset(scaled, idents = idents)
    levels(ThreeCellTypes) <- sort(levels(ThreeCellTypes))
    DoHeatmap(ThreeCellTypes, features=top_n, draw.lines = TRUE)
  }, height = 1000, width = 900)
  
  output$DotPlot_m <- renderPlot({
    input$refresh_4
    # data is what reduction you are using
    types <- c()
    if (isolate(input$Basal_m)){types <- append(types, "Basal")}
    if (isolate(input$Ciliated_m)){types <- append(types, "Ciliated")}
    if (isolate(input$Intermediate_m)){types <- append(types, "Intermediate")}
    if (isolate(input$Preciliated_m)){types <- append(types, "Preciliated")}
    if (isolate(input$Proliferating_m)){types <- append(types, "Proliferating_Basal")}
    if (isolate(input$Secretory_m)){types <- append(types, "Secretory")}
    validate(
      need(length(types) != 0, "Must have at least one cell type")
    )
    treatment <- isolate(input$treatment_m)
    n <- isolate(input$n_m)
    include_downregulated <- isolate(input$include_downregulated_m)
    include_upregulated <- isolate(input$include_upregulated_m)
    validate(
      need((include_downregulated != TRUE || include_upregulated != TRUE), "Cannot plot only downregulated and only upregulated genes at the same time")
    )
    subsetted <- data;
    Idents(subsetted)<-subsetted@meta.data$Treatment
    subsetted <- subset(subsetted, idents = treatment)
    Idents(subsetted)<-subsetted@meta.data$cellType
    markers <- FindMarkers(subsetted, ident.1 = types)
    if (include_downregulated == TRUE){
      markers <- subset(markers, (markers$avg_log2FC < 0))
    }else if (include_upregulated == TRUE){
      markers <- subset(markers, (markers$avg_log2FC > 0))
    }
    markers <- row.names(markers)[1: n]
    validate(
      need(markers[1] != "NA", "It looks like no differentially expressed genes were found")
    )
    # draw the FeaturePlot
    DotPlot(subsetted, features = markers) + coord_flip() + ggtitle(treatment)
  }, height = 1000, width = 900)
  
  output$docs <- renderImage({
    list(src = "../MarKDOWn.png", contentType = "image/png", width="150%")
  }, deleteFile = FALSE)
  
  output$downloadPlot1 <- downloadHandler( filename = function() { paste("plotOutput", '.png', sep='') },
                                           content = function(file) {ggsave(file, plot = last_plot(), device = "png", height = 10, width = 12, limitsize = FALSE)})
  output$downloadPlot2 <- downloadHandler( filename = function() { paste("plotOutput", '.png', sep='') },
                                           content = function(file) {ggsave(file, plot = last_plot(), device = "png", height = plotHeight()/75, width = 12, limitsize = FALSE)})
  output$downloadPlot3 <- downloadHandler( filename = function() { paste("plotOutput", '.png', sep='') },
                                           content = function(file) {ggsave(file, plot = last_plot(), device = "png", height = plotHeight2()/75, width = 12, limitsize = FALSE)})
  output$downloadPlot4 <- downloadHandler( filename = function() { paste("plotOutput", '.png', sep='') },
                                           content = function(file) {ggsave(file, plot = last_plot(), device = "png", height = plotHeight2()/75, width = 12, limitsize = FALSE)})
  output$downloadPlot5 <- downloadHandler( filename = function() { paste("plotOutput", '.png', sep='') },
                                           content = function(file) {ggsave(file, plot = last_plot(), device = "png", height = plotHeight()/75, width = 12, limitsize = FALSE)})
  output$downloadPlot6 <- downloadHandler( filename = function() { paste("plotOutput", '.png', sep='') },
                                           content = function(file) {ggsave(file, plot = last_plot(), device = "png", height = 14, width = 12, limitsize = FALSE)})
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("plotData", ".csv", sep = "")
    },
    content = function(file) {
      if (!input$view_raw_h){}
      write.csv(ggplot_build(last_plot())$plot$data, file, row.names = FALSE)
    }
  )
  
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste("plotData", ".csv", sep = "")
    },
    content = function(file) {
      if (!input$view_raw_m){}
      write.csv(ggplot_build(last_plot())$plot$data, file, row.names = FALSE)
    }
  )
  
  output$downloadData3 <- downloadHandler(
    filename = function() {
      paste("plotData", ".csv", sep = "")
    },
    content = function(file) {
      if (!input$view_raw_d){}
      write.csv(ggplot_build(last_plot())$plot$data, file, row.names = FALSE)
    }
  )
  
  output$table = DT::renderDataTable({
    if (!input$view_raw_h){}
    p <- last_plot()
    ggplot_build(p)$plot$data
  })
  output$table2 = DT::renderDataTable({
    if (!input$view_raw_d){}
    p <- last_plot()
    ggplot_build(p)$plot$data
  })
  output$table3 = DT::renderDataTable({
    if (!input$view_raw_m){}
    p <- last_plot()
    ggplot_build(p)$plot$data
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
