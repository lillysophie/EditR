# EditR UI
editrversion <- "1.0.10" 

shinyUI(
  
  fluidPage(

  titlePanel(list(HTML('<a href="https://twin-cities.umn.edu/"><img src=https://upload.wikimedia.org/wikipedia/commons/thumb/f/f9/Minnesota_Golden_Gophers_logo.svg/2560px-Minnesota_Golden_Gophers_logo.svg.png class="img.unframed" align="right" alt="UMN" height="95" width="170"></a>'), 
                  paste("EditR", editrversion)),
             windowTitle="EditR: Edit Deconvolution by Inference of Traces in R"),
  
  sidebarPanel(
    tags$head( tags$link(rel="stylesheet", type="text/css", href="styleupdates.css")),
    tags$img(href = "https://www.cancer.umn.edu/bio/osteosarcoma-staff/branden-moriarity", src="https://image.ibb.co/bStZTb/editr_logo_v3.png", alt="EditR", width="90%", class="unframed", align="center"),
    checkboxInput("example", "Load Example Data", FALSE),
    conditionalPanel(condition = "!input.example",
                     p("Please upload your Sanger sequencing file(s) here"),
                     fileInput(inputId = "file", 
                               label = 'Upload .ab1 File of edited sample(s)',
                               multiple = TRUE),
                     fileInput(inputId = "file_nt",
                               label = 'Upload .ab1 File of non-treated sample',
                               multiple = FALSE),
                     p("Please enter in your ~20bp guide sequence in 5'-3' orientation"),
                     textInput(inputId = 'guide',
                               label = 'Enter gRNA sequence'),
                     checkboxInput("guide.is.reverseComplement",
                                   "Guide sequence is reverse complement", FALSE),
    
                     p("Manually enter the 5' start and 3' end of the good sequencing below to trim the sanger sequening read. To determine which values to use, click the 'Data QC' tab and cursor over the 'Data QA: Signal and noise plot'."),
                     p("For example, if clean sequencing starts at the 56th bp of the chromatogram and ends at 405th bp, enter those respective values."),
                     numericInput(inputId ="trim5", label = "5' start", value = NA),
                     numericInput(inputId ="trim3", label = "3' end", value = NA),
                     p("Note: aggressive trimming may lead to inaccurate P-values, make sure to trim to the region that represents the majority of your sample."),
                     numericInput("pvalcutoff", "P-value cutoff:", 0.01, min = 0, max = 1),
                     p("Please enter your cutoff probability (p-value for significance), see instructions for details"),
                     
                     checkboxInput("heterozygous", "Heterozygous mutation", FALSE),
                     conditionalPanel(condition = "input.heterozygous",
                                      numericInput(inputId = "het_position",
                                                   label = "Heterozygous position in guide (1 = first base)",
                                                   value = NA, min = 1, step = 1),
                                      p("At this position the Abbott formula, (edited - non-treated) / (100 - non-treated) * 100, is used to calculate percent editing. It requires a non-treated sample.")
                     )
        ),
    p("Optionally enter the position within the guide (1 = first base of the guide) where editing is expected. It will be highlighted in the chromatogram."),
    numericInput(inputId = "expected.edit.position",
                 label = "Expected editing position in guide",
                 value = NA, min = 1, step = 1),
      p("Data are not stored on the server after the analysis is ran.")
  ),
  
  mainPanel(
    
    tabsetPanel(type = "tabs",
                
                #### Instructions tab
                #source("./instructions-content.R", local = TRUE)$value,
                # tabPanel("Instructions",
                # includeHTML("instructions.html")),
                
                tabPanel("Instructions",
                         includeMarkdown("instructions.md")),
                
                #### Data QC tab
                tabPanel("Data QC",
                         p(strong("Data QC is shown for every uploaded edited sample."), " Each sample has four plots:"),
                         tags$ul(
                           tags$li(strong("Total peak area before filtering:"), " the total peak area (A + C + G + T) at each base position before QC filtering. The shaded region denotes the guide RNA."),
                           tags$li(strong("Signal and noise plot:"), " the peak area of the filtered dataset, showing regions with high noise that might affect the sensitivity of the predicted editing. The first 20 bases are removed as they are generally poor quality, then positions whose total peak area is below one tenth of the mean are removed. This mainly trims the low quality tail of the sequencing run (unless trim values are entered)."),
                           tags$li(strong("Percent noise peak area:"), " the noise at each base position as the percent of the total peak area. The shaded region is where the gRNA matches, which might have a peak due to base editing."),
                           tags$li(strong("Guide RNA region chromatogram:"), " the chromatogram of the guide region. The shaded region shows where a double peak was detected during basecalling, although editing can still occur if a double peak was not detected.")
                         ),
                         textOutput(outputId = "trimmedrange"),
                         hr(),
                         uiOutput("qc_sections")
                ),
                
                #### Predicted editing tab
                tabPanel("Predicted Editing",
                         h2("Edited sample(s)"),
                         uiOutput("editing_summary"),
                         uiOutput("edited_sections"),
                         p("The top plot is the chromatogram of the protospacer. The boxed column marks the expected editing position. The bottom plot shows the percent area of the signal for each base (ACGT) at each position along the guide. Bases that are significantly different from the noise are colored in, color coded relative to their percent area. Most positions of the guide only have one base colored in, as this is due to only one peak being present at that position. In the case of editing, there are more than one base colored in at a position."),
                         conditionalPanel(condition = "!input.example",
                           h2("Non-treated control"),
                           uiOutput("nt_section")
                         ),
                         conditionalPanel(condition = "input.heterozygous && !input.example",
                           h2("Normalised Base Areas"),
                           uiOutput("norm_sections"),
                           p("Base areas of each edited sample normalised to the non-treated sample. At the heterozygous position the Abbott formula is used; at all other positions the non-treated signal is subtracted (provisional). The chromatogram is an idealised reconstruction with one peak per base and position.")
                         ),
                         h2("Base info"),
                         p("You can use the information in the following table to better understand the values given for the percent
                           area in the 'Table plot'. The average percent signal shows you the average percent area (area of signal / all
                           other base areas), and (100 - average percent signal) shows you how much noise is present during the reading 
                           of that base. The amount of noise present in each base is also reflected in the 'Critical percent value', 
                           which is the critical value for the percent area for that base -- any percent area value above that would 
                           be called significantly different from the noise, and therefore editing has occurred. This critical value is 
                           calculated based on the P-value cutoff specified (default is 0.01). The model mu is the mu term for the zero
                            adjusted gamma distribution to model the noise -- a higher value of mu means that the predicted editing value
                            is less accurate. Filliben's correlation is how well the noise is modelled by a zero adjusted gamma 
                           distribution. Lower values of Filliben's correlation means less confidence in EditR for predicting editing,
                           your values should be above 0.90, if not, this means that your sequencing file likely has problems
                           with data quality."),
                         uiOutput("baseinfo_sections"),
                         h2("Quad plot"),
                         uiOutput("quad_sections"),
                         p("This plot shows the percent area of the background bases (the bases that 
                           were not the guide squence), the line denotes the crtical value cutoff based on the P-value
                           cutoff specified (default is 0.010.")
                ),                
                
                
                #### Report download tab
                tabPanel("Download Report",
                         h2("Download PDF report here"),
                         p("Click the button below to download a PDF report with the results for all uploaded samples (every edited sample and the non-treated sample). For each sample the report contains the Data QC plots, the guide chromatogram with the percent area table, the base info table and the quad plot. If the heterozygous mutation option is ticked, the normalised base areas of each edited sample are included as well. Please be patient, as the report may take some time to generate."),
                         downloadButton('downloadReport')
                         )
                
                )
    

)
  
  )
)
