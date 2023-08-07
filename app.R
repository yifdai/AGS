library(shiny)
library(Hmisc)
library(DT)
library(png)
library(shinythemes)
library(shinyjs)
library(data.table)
require(rrBLUP)
library(dplyr)
library(pryr)

options(shiny.maxRequestSize = 1024*1024^2)

setwd("/opt/app/Geneselect")  # Change to user dir in the docker image.

load(file = "./example/rrblup/data/example.RData")  # Change to user dir in the docker image.

ui <- fluidPage(

  titlePanel(
    tags$div(style = "margin-top: -22px;",  
             tags$img(src = "title1.png", height = "50px", width = "400px"),
             ""
    )
  ),
  theme = shinytheme("readable"),
  
  tabPanel("Home",
           fluidPage(
             useShinyjs(),
             sidebarLayout(
               sidebarPanel(
                 selectInput(inputId = "model", label = "Select model",
                             choices = c("rrblup", "DNNGP3", "BA", "BB", "BC", "BL","BRNN","BRR","EGBLUP","RF","EN", "GBLUP", "LASSO", "RKHS", "RR", "SVM")),
                 
                 checkboxInput(inputId = "useExample", label = "Use example data", value = F),
                 
                 # selectInput(inputId = "filetypeselect", label = "Select what type of raw data you have",
                 #             choices = c("Genotype and phenotype data, kinship data all stored in one RData file", "Genotype and phenotype data stored in csv file, but kinship data stored in RData file.")),
                 

                 fileInput(
                 inputId = "genodata",
                 label = tags$div(
                   "Upload Genotype Data File",
                   tags$br(),
                   tags$span(style = "color: grey; font-size: 0.9em;", "Here is for Genotype data input.")
                 ),
                 accept = c(".vcf", ".vcf.gz")),
                 
                 fileInput(
                   inputId = "phenodata",
                   label = tags$div(
                     "Upload Phenotype Data File",
                     tags$br(),
                     tags$span(style = "color: grey; font-size: 0.9em;", "Here is for Phenotype data input.")
                   ),
                   accept = c(".txt", ".csv")),

                 
                 conditionalPanel(condition = "input.useExample == true && input.model == 'rrblup'",
                                  selectInput(inputId = "rrbluptrait", label = "Select rrblup example trait",
                                              choices = c("TC","SN","RS","TS","PC","LL","LN","LW","Nic","TN","IC","FGT","RBW","RBSH"))
                 ),
                 
                 conditionalPanel(condition = "input.useExample == false",
				                          conditionalPanel(condition = "input.model == 'DNNGP3'",
				                                           numericInput(inputId = "batch_size", tags$abbr("Batch size", title = "Training batch size"), value = 32),
				                                           numericInput(inputId = "lr", tags$abbr("Learning rate", title = "Initial learning rate"), value = 0.001),
				                                           numericInput(inputId = "epoch", tags$abbr("Number of epochs", title = "Number of training iterations"), value = 100),
				                                           numericInput(inputId = "dropout1", tags$abbr("First Dropout", title = "Dropout rate after the first layer"), value = 0.2),
				                                           numericInput(inputId = "dropout2", tags$abbr("Second Dropout", title = "Dropout rate after the second layer"), value = 0.2),
				                                           numericInput(inputId = "patience", tags$abbr("Patience for learning rate reduction", title = "Number of epochs with no improvement after which learning rate will be reduced"), value = 10),
				                                           numericInput(inputId = "seed", tags$abbr("Random Seed", title = "Seed for random number generator"), value = 123),
				                                           numericInput(inputId = "cv", tags$abbr("Number of folds for cross-validation", title = "Number of cross-validation folds"), value = 5),
				                                           numericInput(inputId = "part", tags$abbr("Part for validation set", title = "Fold to be used as the validation set"), value = 1),
				                                           numericInput(inputId = "earlystopping", tags$abbr("Early stopping threshold", title = "Number of epochs with no improvement after which training will be stopped"), value = 10),
				                                           numericInput(inputId = "pcanum", tags$abbr("PCA number", title = "The PCA number for converting the file to .tsv."), value = 50)
				                          ),
				                          conditionalPanel(condition = "input.model != 'DNNGP3' && input.model != 'rrblup' && input.model != 'DeepGS'",
				                                           selectInput("impute_method", "Impute Method", choices = c("mni", "emi"), selected = "mni"),
				                                           numericInput("MAXNA", "Max NA", value = 0.2, min = 0, max = 1, step = 0.01),
				                                           numericInput("MAF", "MAF", value = 0.05, min = 0, max = 1, step = 0.01),
				                                           textInput("reduct.size", "Reduct Size", "NULL"),
				                                           textInput("r2", "R2", "NULL"),
				                                           textInput("pval", "P Value", "NULL"),
				                                           textInput("MAP", "MAP", "NULL"),
				                                           checkboxInput(inputId = "trimMarkers", label = "Trim markers", value = F),
									   conditionalPanel(condition = "input.trimMarkers == true",
                                                                                            numericInput(inputId = "trimnum", tags$abbr("Number of Markers to involve", title = "Number of Markers you want to keep"), value = 100000)
                  ),
				                          )

                 ),
                  
		  actionButton(inputId = "GMI", label = "Generate match ids phenotype file"),	    
	          actionButton(inputId = "runAnalysis", label = "Run Analysis")
                 
               ),
               
      				 mainPanel(
      				   textOutput("runStatus"),
				         uiOutput("dynamic_trait_select"),
      				   verbatimTextOutput("txtOut"),
      				   actionButton("show_plot_table", "Show Plot and Table"),
      				   uiOutput("dynamic_output")
               )
             )
           )
  )
  
  # tabPanel("Help",
  #          fluidPage(
  #            h1("Help Page"),
  #            p("This is a help page.")
  #          )
  # )
  
)

server <- function(input, output, session) {
  
  values <- reactiveValues(bwgs_result = NULL,
                           rrblup_result = NULL,
                           dnngp_result = NULL)
  
  observeEvent(input$runAnalysis, {
    values$bwgs_result <<- FALSE
    values$rrblup_result <<- FALSE
    values$dnngp_result <<- FALSE
  })

  source("/opt/app/Geneselect/model/plot.R")

  runStatus <- reactiveVal("This is a R shiny app for Gene selection. Here is for the basic instructions. The running process would take 1-10 minutes depends on your data file size. After input the correct input files, the butten 'Run Analysis' will pop up.")
  output$runStatus <- renderText({
    runStatus()
  })
  
  shinyjs::disable("show_plot_table")
  shinyjs::disable("runAnalysis")
  shinyjs::disable("GMI")
  
  selected_trait <- reactiveVal()
  
  observeEvent(input$genodata, {
    ext_genotype <<- tools::file_ext(input$genodata$name)
    if (ext_genotype == "vcf" || ext_genotype == "gz" && grepl("\\.vcf\\.gz$", input$genodata$name)) {
      showNotification(paste0("Upload successful. The input of the Genotype data file is a VCF file. The sample name is ", input$genodata$name))
      # Enable the button if the file is correct
    } else {
      showNotification("Upload failed. The input of the Genotype data file is not a VCF file, please try again.")
      # Disable the button if the file is incorrect
    }
  })
  
  observeEvent(input$phenodata, {
    ext_phenotype <<- tools::file_ext(input$phenodata$name)
    if (ext_phenotype == "txt") {
      showNotification(paste0("Upload successful. The input of the Phenotype data file is a txt file. The sample name is ", input$phenodata$name))
      pheno_dt <<- fread(input$phenodata$datapath, na.strings = c("NA", "na", "NaN", "NONE", "None", "none", "null", "Null", "NULL", ""))
      pheno_dt <<- as.data.frame(pheno_dt)
      cat(paste0(sum(is.na(pheno_dt)), "\n"))
      file_type <<- " "
      shinyjs::enable("runAnalysis")
      shinyjs::enable("GMI")
      # Enable the button if the file is correct
    } else if (ext_phenotype == "csv"){
      showNotification(paste0("Upload successful. The input of the Phenotype data file is a csv file. The sample name is ", input$phenodata$name))
      pheno_dt <<- fread(input$phenodata$datapath, na.strings = c("NA", "na", "NaN", "NONE", "None", "none", "null", "Null", "NULL", ""))
      pheno_dt <<- as.data.frame(pheno_dt)
      file_type <<- ","
      shinyjs::enable("runAnalysis")
      shinyjs::enable("GMI")
      # Enable the button if the file is correct
    }else {
      showNotification("Upload failed. Input of the Phenotype data file is not a txt file, please try again.")
      # Disable the button if the file is incorrect
    }

    pheno_traits <- colnames(pheno_dt)[-1]
    
    # Render the selectInput with the traits
    output$dynamic_trait_select <- renderUI({
        selectInput("selected_trait", "Select a trait:", choices = pheno_traits)
    })
  })
  
  observeEvent(input$GMI, {
    genofilename_no_ext <<- sub("\\.vcf(\\.gz)?$", "", input$genodata$name)
    phenofilename_no_ext <<- gsub("\\.txt|\\.csv", "", input$phenodata$name)
    system(paste0("bash /opt/app/Geneselect/opt/convert.sh ", genofilename_no_ext, " ", input$genodata$datapath))
    genotxt_file_dir <<- paste0("/mount/result/", genofilename_no_ext, "/")
    vcf <- fread(paste0(genotxt_file_dir, genofilename_no_ext, ".fam"), header = FALSE, sep = " ")
    sid <- colnames(pheno_dt)[1]
    
    pheno <- pheno_dt
    
    common_names <- intersect(pheno[[sid]], vcf[[2]])
    
    if (length(common_names) > 0) {
      cat("Original ids in pheno: \n")
      cat(pheno[[sid]], "\n")
      cat("Original ids in vcf: \n")
      cat(vcf[[2]], "\n")
      
      cat("Matching names in vcf iid: \n")
      cat(common_names, "\n")
      
      pheno_new <- pheno %>%
        filter(pheno[[sid]] %in% vcf[[2]]) 
      
      # Now reorder pheno_new based on the order in vcf[[2]]
      pheno_new <- pheno_new[match(vcf[[2]], pheno_new[[sid]]), ]
      
      pheno_new <- pheno_new %>% filter(!is.na(pheno_new[[sid]]))
      
      write.csv(pheno_new, paste0(genotxt_file_dir, phenofilename_no_ext, "_matched.csv"), row.names = F)
      cat("--------------------------Get the matched pheno data----------------------------\n")
      cat("You can now check the matched file in /result folder.\n")
    } else {
      cat("Error: No common individuals in vcf file and phenotype file, please check your input.\n")
    }
  })

  observeEvent(input$selected_trait, {
    selected_trait(input$selected_trait)
    cat(input$selected_trait)
  })

  output$dynamic_output <- renderUI({
    if (input$model != "DNNGP3" && input$model != "rrblup" && input$model != "DeepGS") {  
      tabsetPanel(
        tabPanel("Plot", 
                 fluidRow(
                   column(6, plotOutput("bwgs_plot1")),
                 )
        ),
        tabPanel("Table",
                 tableOutput("bwgs_table1"),
                 tableOutput("bwgs_table2")
        )
      )
    } else if (input$model == "rrblup") {  
      tabsetPanel(
        tabPanel("Plot", 
                 fluidRow(
                   column(6, plotOutput("rrblup_plot1")),
                 ),
                 fluidRow(
                   column(6, plotOutput("rrblup_plot2"))
                   # column(6, plotOutput("rrblup_plot3"))
                 )
        ),
        tabPanel("Table",
                 tableOutput("rrblup_table"),
		 tableOutput("rrblup_table2")
        )
      )
    } else if (input$model == "DNNGP3") {
      tabsetPanel(
        tabPanel("Plot", 
                 fluidRow(
                   column(6, plotOutput("dnngp_plot1")),
                 ),
                 fluidRow(
                   column(6, plotOutput("dnngp_plot2")),
                 )
        ),
        tabPanel("Table",
                 tableOutput("dnngp_table1")
        )
      )
    }
  })
  
  observeEvent(input$show_plot_table, {
    
    if (values$rrblup_result) {
      cat("rrblup method done\n")
      cat(paste0("Load path: ", genotxt_file_dir_rrblup, "\n"))
      output$rrblup_plot1 <- renderImage({
        list(src = paste0(genotxt_file_dir_rrblup,"traits-", traits_suf,".png"), contentType = "image/png", height = 400)
      }, deleteFile = FALSE)
      
      # output$rrblup_plot2 <- renderImage({
       # req(selected_trait()) # Ensure the selected trait is available
      #  list(src = paste0(genotxt_file_dir_rrblup, "plot_", selected_trait(), ".png"), contentType = "image/png", height = 400)
     # }, deleteFile = FALSE)
      
      output$rrblup_plot2 <- renderImage({
        req(selected_trait()) # Ensure the selected trait is available
        list(src = paste0(genotxt_file_dir_rrblup, selected_trait(),"/plot2_", selected_trait(), ".png"), contentType = "image/png", height = 400)
      }, deleteFile = FALSE)
      
      output$rrblup_table <- renderTable({
        output2
      })

      output$rrblup_table2 <- renderTable({
	output_pred_rrblup <- read.csv(paste0(genotxt_file_dir_rrblup, selected_trait(), "/Prediction.csv"))
	output_pred_rrblup
      })
      
    } else if (values$dnngp_result) {
      cat("dnngp method done\n")
      cat(paste0("Load path: ", dnngp_output, "\n"))
      output$dnngp_table1 <- renderTable({
        va <- read.csv(paste0(dnngp_output, selected_trait(), "/Prediction.ALL.csv"))
        va
      })
      
      output$dnngp_plot1 <- renderImage({
        req(selected_trait()) # Ensure the selected trait is available
        list(src = paste0(dnngp_output, selected_trait(),"/loss_over_time_", selected_trait(),".png"), contentType = "image/png", height = 400)
      }, deleteFile = FALSE)
      
      output$dnngp_plot2 <- renderImage({
        req(selected_trait()) # Ensure the selected trait is available
        list(src = paste0(dnngp_output, selected_trait(), "/validation_plot_", selected_trait(),".png"), contentType = "image/png", height = 400)
      }, deleteFile = FALSE)
      
    } else if (values$bwgs_result) {
      cat(paste0(input$model ," method done\n"))
      cat(paste0("Load path: ", bwgs_output, "\n"))
      
      bwgs_train <- read.csv(paste0(bwgs_output, selected_trait(), "_", input$model,"_summary_cv.csv"))
      bwgs_prediction <- read.csv(paste0(bwgs_output, selected_trait(), "_", input$model,"_prediction.csv"))
      
      output$bwgs_table1 <- renderTable({
        bwgs_train  # train_summary_table
      })
      output$bwgs_table2 <- renderTable({
        bwgs_prediction  # prediction_table
      })
      
      output$bwgs_plot1 <- renderImage({
        req(selected_trait()) # Ensure the selected trait is available
        list(src = paste0(bwgs_output, "validation_plot_", input$model, "_", selected_trait(),".png"), contentType = "image/png", height = 400)
      }, deleteFile = FALSE)
    }
  })

  #########################################################################
  # rrblup
  observeEvent(input$runAnalysis, {
      req(input$useExample == FALSE)
      if (input$model == "rrblup") {
      	withProgress(message = 'rrblup running', value = 0, {      
            incProgress(1/4, detail = "Data pre-processing")
            
	          sample_ID <- colnames(pheno_dt)[1]
      	    traits <- colnames(pheno_dt)[-1]
      	    acc_result <- data.frame(matrix(ncol = 8, nrow = 0))
      	    
      	    # Obtain the filename without extension / path specification
      	    genofilename_no_ext <<- sub("\\.vcf(\\.gz)?$", "", input$genodata$name)
      	    genotxt_file_dir <<- paste0("/mount/result/", genofilename_no_ext, "/")
      	    genotxt_file_dir_rrblup <<- paste0("/mount/result/", genofilename_no_ext, "/rrblup-result/")
      	    
      	    system(paste0("bash /opt/app/Geneselect/opt/convert.sh ", genofilename_no_ext, " ", input$genodata$datapath))
            	    
      	    cat("------------------------------checking data---------------------------------\n")
      	    sid <- colnames(pheno_dt)[1]
      	    
            # Start a new session for the convert.sh script
            system(paste0("sh /opt/app/Geneselect/opt/convert_rrblup.sh ",genofilename_no_ext, " ",input$genodata$datapath))
            incProgress(1/4, detail = "Data pre-processing step complete")
            cat("------------------------------data preprocessing complete----------------------------\n")          
            source("./model/rrblupf.R")
            source("./model/grmf.R")

      	    cat(paste0("---------------------------trying to get kinship matrix------------------------------\n")) 
            grm_data <- ReadGRMBin(prefix = paste0(genotxt_file_dir, genofilename_no_ext))
            kinship_matrix <- construct_kinship_matrix(grm_data)
            write.csv(kinship_matrix, file = paste0(genotxt_file_dir_rrblup, "/kinship_matrix.csv"), row.names = FALSE)
	          cat(paste0("------------------------------------done---------------------------------------------\n"))
            out <- data.frame("traits"=traits,"r"=1,"sd"=1)
            incProgress(1/4, detail = "Get raw data.")
        
            withProgress(message = 'Running rrblup with user data', value = 0, {
              require(rrBLUP)
              require(Hmisc)
              traits_text <- ""
              Rprof(filename = paste0(genotxt_file_dir, "memory_profile.out"), memory.profiling = TRUE)
              summary_texts <- list()
              for (i in 1:length(traits)) {
                  setProgress(message = paste0('Processing file ', traits[i]))
                  start_time <- Sys.time()
                  
                  trait <- traits[i]
                  summary_texts <- append(summary_texts, paste0("For trait ", trait, ", in total ", sum(is.na(pheno_dt[, trait])), " samples were predicted."))
            		  cat(paste0("---------------------------Processing for trait ", trait, "------------------------------\n"))
            		  traits_text <- paste0(traits_text, substr(trait, 1, 1))
            		  cat(paste0("---------------------------Preparing for Input------------------------------\n"))
                  input_now <- prep_input(pheno_dt, kinship_matrix, colnames(pheno_dt)[1], trait = trait)
                  cat(paste0("---------------------------------done------------------------------\n"))
            		  cat(paste0("---------------------------------CV rrblup----------------------------------\n"))
            		  out_now <- cv_rrblup(y = input_now$y, Z = input_now$Z, t_now=trait, genotxt_file_dir_rrblup)
            		  cat(paste0("---------------------------------done------------------------------\n"))
            		  end_time <- Sys.time()
                  Rprof(NULL)
                  cat(end_time - start_time)
                  mem_summary <- summaryRprof(paste0(genotxt_file_dir, "memory_profile.out"))

        	        r <- mean(out_now,na.rm = T)
                  sd <- sd(out_now,na.rm = T)
                  out[i,2:3] <- c(r,sd)
                  
                  system(paste0("mkdir -p ", genotxt_file_dir_rrblup, trait))

                  load(paste0(genotxt_file_dir_rrblup, trait, ".RData"))
            		  cat(paste0("--------------------------------writing results-------------------------------\n"))
                  # Get predicted values for missing y
            		  rr <- cor(pre1, pre2)
                  # create_plot(Z %*% est$u, y, trait, paste0(genotxt_file_dir_rrblup, "plot_", trait, ".png"))
                  create_plot(pre1, pre2, trait, paste0(genotxt_file_dir_rrblup, trait,"/plot_validation_", trait, ".png"))
                  
                  incProgress(1/length(traits)) # Update progress bar
              }
            
              incProgress(1/4, detail = "rrblup done, generating plots.")
              require(Hmisc)
              # Calculate standard deviation
              out_sub <- out
              out_sub$sd <- out_sub$r*0.02
              
              # Create dataframe for ggplot
              data <- data.frame(
                x = 1:nrow(out_sub),
                y = out_sub$r,
                y_min = out_sub$r - out_sub$sd,
                y_max = out_sub$r + out_sub$sd,
                label = as.character(out_sub$traits)
              )
              
              # Create plot
              p <- ggplot(data, aes(x = x, y = y)) +
                geom_point(color="#1f77b4", size=1.5, alpha=0.6) +
                geom_errorbar(aes(ymin=y_min, ymax=y_max), width=0.2, color="blue") +
                geom_hline(aes(yintercept=mean(y)), linetype="dashed", color="red") +
                labs(x="", y="Accuracy", title="traits") +  # remove caption from labs
                theme_minimal() +
                theme(plot.title = element_text(hjust = 0.5), 
                      text = element_text(size=14),
                      axis.title = element_text(face="bold"),
                      legend.position = "none") +
                scale_x_continuous(breaks = data$x, labels = data$label)
              
              # Save plot
              ggsave(filename = paste0(genotxt_file_dir_rrblup,"traits-", traits_text,".png"), plot=p, width=5, height=5, dpi=300)

              out <- data.frame("traits"=traits,"r"=1,"sd"=1)
              output2 <<- loop_cv(pheno_dt, kinship_matrix, colnames(pheno_dt)[1], out, acc_result, genotxt_file_dir_rrblup)
              cat(paste0("-------------------------------get the acc.table------------------------------\n"))
	      write.table(output2, file=paste0(genotxt_file_dir_rrblup,"/acc.table"),row.names = F)
	    })
      traits_suf <<- traits_text
  
	    plots <- c(paste0(genotxt_file_dir_rrblup, "traits-", traits_text,".png"), 
	               # paste0(genotxt_file_dir_rrblup, selected_trait(),"/plot_", selected_trait(), ".png"), 
	               paste0(genotxt_file_dir_rrblup, selected_trait(),"/plot2_", selected_trait(), ".png"))
	    plots_exists <- all(file.exists(plots))
	    # cat(file.exists(plots))

	    # Debugging step 1
	    print(dim(summary_texts))
	    print(summary_texts)
	    str(summary_texts) # Check the structure of summary_texts

	    # Debugging step 2
	    for(i in 1:length(summary_texts)){
		    print(summary_texts[[i]])  # print each element separately
	    }

      	    for(item in unlist(summary_texts)) {
      	      cat(item, "\n")  
      	      cat("-----------------------------------------------\n") 
      	    }
      	    
      	    if (plots_exists) {
      	      values$rrblup_result <<- TRUE
      	      shinyjs::enable("show_plot_table")
      	      cat("RRBLUP_DONE\n")
	      cat("You can check the result on the server or check the result under /rrblup_result folder\n")
	            runStatus(paste0("rrblup completed, the result can be found in ", genotxt_file_dir_rrblup))
      	    } else {
      	      cat("RRBLUP_NULL\n")
	      warning("Error occured, the plots are not correctly generated")
      	    }
      	})
        }
  })
  ###########################################################
  
  ###########################################################
  # DeepGS
  observeEvent(input$runAnalysis,{
    req(input$useExample == FALSE)
    if (input$model == "DeepGS") {
      withProgress(message = "DeepGS running", value = 0, {
      	incProgress(1/3, detail = "Preparing data for DeepGS method")
      	# Obtain the filename without extension
      	genofilename_no_ext <<- sub("\\.vcf(\\.gz)?$", "", input$genodata$name)
        cat(genofilename_no_ext)
      	# Specify the directory of the result
      	genotxt_file_dir <<- paste0("/mount/result/", genofilename_no_ext, "/")
      	system(paste0("mkdir ", genotxt_file_dir))
      	genotxt_file_deepgs <<- paste0("/mount/result/", genofilename_no_ext, "/DeepGS_result")
      	system(paste0("mkdir ", genotxt_file_deepgs))
      	cat(genotxt_file_deepgs)
      	dp_command <- paste0("Rscript /opt/app/Geneselect/model/DeepGS_get_raw.R ", input$genodata$datapath, " ", input$phenodata$datapath, " ", genotxt_file_deepgs, "/", genofilename_no_ext)
      	cat(dp_command)
        system(dp_command)
        incProgress(1/3, detail = "Running DeepGS method")
      	# Wait for DeepGS to finish
      	while (TRUE) {
      		log_lines <- readLines("/mount/mnt.log")
      		if ("DeepGS finished" %in% log_lines) {
      			break
      		}
      		Sys.sleep(1)
      	}
      	
      	incProgress(1/3, detail = "DeepGS finished")
      	file.remove("/mount/mnt.log")
     })
    }
  }) 
  ########################################################

  ###########################################################
  # DNNGP3
  observeEvent(input$runAnalysis, {
    req(input$useExample == FALSE)
    if (input$model == "DNNGP3") {
      withProgress(message = "DNNGP3 running", value = 0, {
        
        sid <- colnames(pheno_dt)[1]
        incProgress(1/3, detail = paste0("Running plink commands"))  
        # Get raw data info
        genofilename_no_ext <- sub("\\.vcf(\\.gz)?$", "", input$genodata$name)
        genotxt_file_dir <- paste0("/mount/result/", genofilename_no_ext, "/")
        genoeigenvec_path <- paste0(genotxt_file_dir, "pca.eigenvec")
        
        initial_memory <- mem_used()
        
        # Rprof(filename = paste0(genotxt_file_dir, "memory_profile.out"), interval = 0.01, memory.profiling = TRUE)
        start_time <- Sys.time()
        
        system(paste0("bash /opt/app/Geneselect/opt/convert.sh ", genofilename_no_ext, " ", input$genodata$datapath))
        system(paste0("/opt/miniconda3/envs/DNNGP3/bin/python /opt/app/Geneselect/opt/plink2dnngp_pca.py ", genoeigenvec_path, " ", genoeigenvec_path))
        
      	cat(paste0("---------------------------get pca.eigenvec------------------------------\n"))
      	
      	incProgress(1/3, detail = paste0("DNNGP3 processing"))
        
        sid <- colnames(pheno_dt)[1]
        geno_names <- read.table(paste0(genotxt_file_dir, genofilename_no_ext,".fam"), sep = " ", header = FALSE)
        colnames(geno_names) <- c("FID", "IID", "V1", "V2", "V3", "V4")
        pheno_ids <- pheno_dt[, sid]
        geno_ids <- geno_names[, "IID"]
        
        common_ids <- intersect(pheno_ids, geno_ids)
        num_common <- length(common_ids)
        
        num_pheno <- length(pheno_ids)
        num_geno <- length(geno_ids)
        min_num <- min(num_pheno, num_geno)
        
        if(num_common < min_num){
          cat(paste0("In total ", num_common, " common ids.\n"))
        }
        
        # Calculate the length difference between the two id lists
        diff_length <- length(pheno_ids) - length(geno_ids)
        
        # If pheno_ids is shorter, add NAs to make it the same length as geno_ids
        if (diff_length < 0) {
          pheno_ids <- c(pheno_ids, rep(NA, abs(diff_length)))
        }
        
        # If geno_ids is shorter, add NAs to make it the same length as pheno_ids
        if (diff_length > 0) {
          geno_ids <- c(geno_ids, rep(NA, abs(diff_length)))
        }
        
        # Now, both vectors should be of the same length, and we can create the data frame
        df <- data.frame("pheno_ids" = pheno_ids, "geno_ids" = geno_ids)
        compare_ids_path <- paste0(genotxt_file_dir, "compare_ids.csv")
        write.csv(df, file = compare_ids_path, row.names = FALSE)
        cat(paste0("The comparison between pheno ids and geno ids (IID) can be found in ", compare_ids_path, "\n"))
        
        if(length(pheno_ids) != length(geno_ids) || num_common == 0){
          warning("The ids not matching. You should check compare_ids.csv to see if the pheno id and IID in VCF file are corresponding to each other, if not please make them correspond to each other.")
        }
        validate(
          need(length(pheno_ids) == length(geno_ids) && num_common != 0, 
               "The ids not matching. You should check compare_ids.csv to see if the pheno id and IID in VCF file are corresponding to each other, if not please make them correspond to each other.")
        )
        
        summary_texts <- list()
	withProgress(message = "Processing", value = 0, {
	for (trait in colnames(pheno_dt)[-1]) {
	    incProgress(1/length(colnames(pheno_dt)[-1]), detail = paste0("trait ", trait))
            dnngp_out <- paste0(genotxt_file_dir, "DNNGP_result/", trait, "/")
	    
      	    cat(paste0("---------------------------Processing for trait ", trait, "------------------------------\n"))
      	    
      	    system(paste0("mkdir -p ", dnngp_out))
      	    system(paste0("sh /opt/app/Geneselect/opt/convert_dnngp.sh ",genofilename_no_ext, " ",input$genodata$datapath, " ", input$phenodata$datapath," ", sid))

            withProgress(message = "Processing", value = 0, {
            system(paste0("sh /opt/app/Geneselect/opt/pca.sh ", genofilename_no_ext, " ", input$pcanum))
          
            incProgress(1/7, detail = "Get eigen vector")
            system(paste0("/opt/miniconda3/envs/DNNGP3/bin/python /opt/app/Geneselect/opt/plink2dnngp_pca_new.py ", genoeigenvec_path, " ", genoeigenvec_path))
            
            predict_ID <- pheno_dt[is.na(pheno_dt[, trait]), sid] # Get the predict IDs based on NA in trait column
            train_ID <- pheno_dt[!is.na(pheno_dt[, trait]), sid] # Get the training IDs based on not NA in trait column
            # Get the training pheno data
            cat("Prediction IDs\n")
      	    cat(paste0(predict_ID, "\n"))
      	    train_dt <- pheno_dt[!is.na(pheno_dt[, trait]), c(sid, trait)]
            write.table(train_dt, file = paste0(dnngp_out, "train_dt.txt"), sep = "\t", row.names = FALSE)
            
            cat("freading eigenvec file ...\n")
            cat("Depends on your marker number, this step could take a few minutes\n")
            pca_data <- fread(genoeigenvec_path, header = TRUE)
            pca_data <- data.frame(pca_data)
      	    colnames(pca_data)[1] <- "#FID"
      	    print(colnames(pca_data))
            # Split the data into prediction part and training part based on the IDs in the trait column
      	    predict_data <- pca_data[pca_data[, "IID"] %in% predict_ID, ]
      	    train_data <- pca_data[pca_data[, "IID"] %in% train_ID, ]
      	    
      	    if (nrow(predict_data) == 0) {
      	        warning("Their's no NA in this trait, thus no prediction will be conduct.")
      	        summary_texts <- append(summary_texts, paste0("For trait ", trait, ", there is no NA, thus no prediction will be done."))
      	   	next
      	    }
      	    summary_texts <- append(summary_texts, paste0("For trait ", trait, ", in total ", nrow(predict_data), " samples were predicted."))
      
      	    # Save predict and train pca.eigenvec
      	    genoeigenvec_predict_path <- paste0(dnngp_out, "predict_pca.eigenvec")
            genoeigenvec_train_path <- paste0(dnngp_out, "train_pca.eigenvec")
            write.table(predict_data, file = genoeigenvec_predict_path, row.names = FALSE, sep = "\t", quote = FALSE)
            write.table(train_data, file = genoeigenvec_train_path, row.names = FALSE, sep = "\t", quote = FALSE)
            
            genopkl_train <- paste0(dnngp_out, "train.pkl")
            genopkl_predict <- paste0(dnngp_out, "predict.pkl")
            
            incProgress(1/7, detail = "Trim eigen vector file")
            
            system(paste0("/opt/miniconda3/envs/DNNGP3/bin/python /opt/app/Geneselect/model/DNNGP-main/Input_files/tsv2pkl.py ", genoeigenvec_train_path, " ", genopkl_train))
            system(paste0("/opt/miniconda3/envs/DNNGP3/bin/python /opt/app/Geneselect/model/DNNGP-main/Input_files/tsv2pkl.py ", genoeigenvec_predict_path, " ", genopkl_predict))
            
            cat("------------------------------data preprocessing complete----------------------------\n")
            incProgress(1/7, detail = "Data pre-processing step completed.")
            
            incProgress(1/7, detail = "Training model.")
            # Capture the user inputs
            batch_size <- input$batch_size
            lr <- input$lr
            epoch <- input$epoch
            dropout1 <- input$dropout1
            dropout2 <- input$dropout2
            patience <- input$patience
            seed <- input$seed
            cv <- input$cv
            part <- input$part
            earlystopping <- input$earlystopping
            snp <- genopkl_train
            pheno <- paste0(dnngp_out, "train_dt.txt")
            output <- dnngp_out
            
            python_cmd_1 <- paste0("/opt/app/Geneselect/model/DNNGP-main/Scripts/dnngp_runner.py --batch_size ", batch_size,
                                 " --lr ", lr,
                                 " --epoch ", epoch,
                                 " --dropout1 ", dropout1,
                                 " --dropout2 ", dropout2,
                                 " --patience ", patience,
                                 " --seed ", seed,
                                 " --cv ", cv,
                                 " --part ", part,
                                 " --earlystopping ", earlystopping,
                                 " --snp '", snp,
                                 "' --pheno '", pheno,
                                 "' --output '", output, "'")

            cat(paste0("---------------Processing for trait ", trait, ", training model------------------------------\n"))
            train_cmd <- paste0("/opt/miniconda3/envs/DNNGP3/bin/python ", python_cmd_1)
            # Run the model training command
            system(train_cmd)
            
            python_cmd_2 <- paste0("/opt/app/Geneselect/model/DNNGP-main/Scripts/Pre_runner.py --Model ", output,
                                 "/training.model.h5 --SNP ", genopkl_predict,
                                 " --output ", output
            )
            incProgress(1/7, detail = "Using trained model to do the prediction")
            
            # Run the model prediction command
            tryCatch({
	    predict_cmd <- paste0("/opt/miniconda3/envs/DNNGP3/bin/python ", python_cmd_2)
	    }, error=function(e) {
		print(e)
	    })
	    cat(paste0("---------------Processing for trait ", trait, ", predicting------------------------------\n"))
	    system(predict_cmd)
            incProgress(1/7, detail = "Complete.")
            
            runStatus(paste0("DNNGP3 running complete. The result saved in the software_installation_dir/result/", genofilename_no_ext, "/DNNGP_result"))
            
            final_memory <- mem_used()
            memory_change <- final_memory - initial_memory
            
            initial_mem_MB <- round(initial_memory / (1024^2), 2)
            final_mem_MB <- round(final_memory / (1024^2), 2)
            mem_change_MB <- round(memory_change / (1024^2), 2)
            
            cat(paste0("Initial memory used: ", initial_mem_MB, " MB\n"))
            cat(paste0("Final memory used:", final_mem_MB, " MB\n"))
            cat(paste0("emory change during the code execution:", mem_change_MB, " MB\n"))
            
            end_time <- Sys.time()
            cat(paste0("Total running time: ", as.numeric(difftime(end_time, start_time, units = "mins")), " minutes\n"))
            
            # Rprof(NULL)
            # mem_summary <- summaryRprof(paste0(genotxt_file_dir, "memory_profile.out"))
            
            incProgress(1/7, detail = "Visualizing result.")
		      
            # Create the prediction / measure plot
            va <- read.csv(paste0(output, "/Prediction.validation.csv"))
            colnames(va)[2] <- paste0(trait, "_prediction")
            va[, trait] <- pheno_dt[match(va[, sid], pheno_dt[, sid]), trait]
            p1 <- create_plot(va[,paste0(trait, "_prediction")], va[,trait], trait, paste0(output, "/validation_plot_", trait,".png"))
            write.csv(va, file = paste0(output, "/Prediction_validation_", trait,".csv"), row.names = FALSE)
            
            # Create the training history plot
            mh <- read.csv(paste0(output, "/Modelhistory.csv"))
            p2 <- ggplot(mh) +
            geom_line(aes(x=epoch, y=loss, color="Training Loss")) +
            geom_line(aes(x=epoch, y=val_loss, color="Validation Loss")) +
                    labs(x="Epoch", y="Loss",
                    title="Training and Validation Loss Over Time",
                      color="Line") +
                    theme_minimal() +
                    theme(plot.title = element_text(hjust = 0.5),
                    text = element_text(size=12))
         
            # Save the plot
            ggsave(filename=paste0(output, "/loss_over_time_", trait,".png"), plot=p2, width=10, height=7, dpi=300)
            })
        }
        })
        dnngp_output <<- paste0(genotxt_file_dir, "DNNGP_result/")
        
        dnngp_plots <- all(file.exists(c(paste0(output, "/Prediction_validation_", trait,".csv"), paste0(output, "/loss_over_time_", trait,".png"))))
	      
        for(item in unlist(summary_texts)) {
          cat(item, "\n")  
          cat("-----------------------------------------------\n") 
        }
        
        if (dnngp_plots) {
          values$dnngp_result <<- TRUE
          shinyjs::enable("show_plot_table")
          cat("DNNGP_TRUE\n")
          cat("You can check the result on the server or check the result under /DNNGP_result folder\n")
	  runStatus(paste0("DNNGP completed, the result can be found in ", dnngp_output))
	} else {
          cat("DNNGP_NULL\n")
	  warning("Error occured, the plots are not correctly generated")
        }
        
       })
    }
  })
  ########################################################

  ########################################################
  # BWGS
  observeEvent(input$runAnalysis, {
    req(input$useExample == FALSE)
    if (input$model != 'rrblup' && input$model != 'DNNGP3' && input$model != 'DeepGS') {
        withProgress(message = paste0(input$model ," running"), value = 0, {
          # Obtain the filename without extension
          genofilename_no_ext <- sub("\\.vcf(\\.gz)?$", "", input$genodata$name)
          # Specify the directory of the generated txt file
          genotxt_file_dir <- paste0("/mount/result/", genofilename_no_ext, "/")
          # Assume the generated txt file has the same name as the original vcf file, but with .txt extension
          genotxt_file_path <- paste0(genotxt_file_dir,"genotype-final-", genofilename_no_ext, ".txt")
          output <- paste0(genotxt_file_dir, "BWGS_result/")
          
          incProgress(1/3, detail = "Data pre-processing.")
          
          system(paste0("bash /opt/app/Geneselect/opt/convert.sh ", genofilename_no_ext, " ", input$genodata$datapath))
          
          mks <- fread(paste0(genotxt_file_dir, genofilename_no_ext,".bim"), sep = " ")
          nrows_mks <- nrow(mks)
          if (input$trimMarkers == TRUE) {
            system(paste0("bash /opt/app/Geneselect/opt/convert_bwgs_trim.sh ", genofilename_no_ext, " ", input$genodata$datapath, " ", input$phenodata$datapath, " ", input$trimnum))
          } else {
            system(paste0("bash /opt/app/Geneselect/opt/convert_bwgs.sh ", genofilename_no_ext, " ", input$genodata$datapath, " ", input$phenodata$datapath))
          }

          source('./model/BWGSf.R')
          source('./model/BWGS_cv.R')
          
          initial_memory <- mem_used()
          # Rprof(filename = paste0(genotxt_file_dir, "memory_profile.out"), interval = 0.01, memory.profiling = TRUE)
          start_time <- Sys.time()
          incProgress(1/3, detail = "BWGS methods processing")
          
          cat("Start 5 fold CV\n")

          tryCatch({
            bwgs_out <- bwgs_cv(genotxt_file_path, input$phenodata$datapath, input$model, genofilename_no_ext,
                             input$impute_method, 
                             input$MAXNA,
                             input$MAF,
                             input$reduct.size,
                             input$r2,
                             input$pval,
                             input$MAP)
          }, error=function(e) {
            print(e)
          })
          
          # sid <- colnames(pheno_dt)[1]
          # geno_names <- read.table(paste0(genotxt_file_dir, genofilename_no_ext,".fam"), sep = " ", header = FALSE)
          # colnames(geno_names) <- c("FID", "IID", "V1", "V2", "V3", "V4")
          # pheno_ids <- pheno_dt[, sid]
          # geno_ids <- geno_names[, "IID"]
          # 
          # common_ids <- intersect(pheno_ids, geno_ids)
          # num_common <- length(common_ids)
          # 
          # num_pheno <- length(pheno_ids)
          # num_geno <- length(geno_ids)
          # min_num <- min(num_pheno, num_geno)
          # 
          # if(num_common < min_num){
          #   cat(paste0("In total ", num_common, " common ids.\n"))
          # }
          # 
          # # Calculate the length difference between the two id lists
          # diff_length <- length(pheno_ids) - length(geno_ids)
          # 
          # # If pheno_ids is shorter, add NAs to make it the same length as geno_ids
          # if (diff_length < 0) {
          #   pheno_ids <- c(pheno_ids, rep(NA, abs(diff_length)))
          # }
          # 
          # # If geno_ids is shorter, add NAs to make it the same length as pheno_ids
          # if (diff_length > 0) {
          #   geno_ids <- c(geno_ids, rep(NA, abs(diff_length)))
          # }
          # 
          # # Now, both vectors should be of the same length, and we can create the data frame
          # df <- data.frame("pheno_ids" = pheno_ids, "geno_ids" = geno_ids)
          # compare_ids_path <- paste0(genotxt_file_dir, "compare_ids.csv")
          # write.csv(df, file = compare_ids_path, row.names = FALSE)
          # cat(paste0("The comparison between pheno ids and geno ids (IID) can be found in ", compare_ids_path, "\n"))
          # if(length(pheno_ids) != length(geno_ids) || num_common == 0){
          #   warning("The ids not matching. You should check compare_ids.csv to see if the pheno id and IID in VCF file are corresponding to each other, if not please make them correspond to each other.")
          # }
          # validate(
          #   need(length(pheno_ids) == length(geno_ids) && num_common != 0, 
          #        "The ids not matching. You should check compare_ids.csv to see if the pheno id and IID in VCF file are corresponding to each other, if not please make them correspond to each other.")
          # )
          # 
          # tryCatch({
          # bwgs_out <- bwgs(genotxt_file_path, input$phenodata$datapath, input$model, genofilename_no_ext,
          #                     input$impute_method, 
          #                     input$MAXNA,
          #                     input$MAF,
          #                     input$reduct.size,
          #                     input$r2,
          #                     input$pval,
          #                     input$MAP)
          # }, error=function(e) {
          #   print(e)
          # })
          
          incProgress(1/3, detail = "Complete")
          end_time <- Sys.time()
          cat(paste0(end_time - start_time), " hours\n")
          
          final_memory <- mem_used()
          memory_change <- final_memory - initial_memory
          
          initial_mem_MB <- round(initial_memory / (1024^2), 2)
          final_mem_MB <- round(final_memory / (1024^2), 2)
          mem_change_MB <- round(memory_change / (1024^2), 2)
          
          cat(paste0("Initial memory used: ", initial_mem_MB, " MB\n"))
          cat(paste0("Final memory used:", final_mem_MB, " MB\n"))
          cat(paste0("memory change during the code execution:", mem_change_MB, " MB\n"))
          
          end_time <- Sys.time()
          cat(paste0("Total running time: ", as.numeric(difftime(end_time, start_time, units = "mins")), " minutes\n"))
          
#           bv_table <- read.csv(paste0(output, selected_trait(), "_", input$model, "_bv.csv"))
#           
#           for (trait in colnames(pheno_dt)[-1]) {
#             create_plot(bv_table[["pheno"]], bv_table[["bv_predict_mean"]], paste0(trait, "_validation_result"), paste0(output, "/validation_plot_", input$model, "_", trait,".png"))
#           }
#           bwgs_output <<- output
#           
#           for(item in unlist(bwgs_out[[3]])) {
#             cat(item, "\n")  
#             cat("-----------------------------------------------\n") 
#           }
# 
#           if (file.exists(paste0(output, "/validation_plot_", input$model, "_", selected_trait(),".png"))) {
#             values$bwgs_result <<- TRUE
#             shinyjs::enable("show_plot_table")
#             cat("BWGS_TRUE\n")
# 	    cat("You can check the result on the server or check the result under /DNNGP_result folder\n")
# 	    runStatus(paste0(input$model, " completed, the result can be found in ", bwgs_output))
#           } else {
#             cat("BWGS_NULL\n")
# 	    warning("Error occured, the plots are not correctly generated")
#           }
          
        })
    }
  })
  
  #####################################################
  
  
  
  
  

  #####################################################
  ## Example part
  observeEvent(input$runAnalysis, {
    req(input$useExample == TRUE)
    if (input$model == "rrblup") {
      runStatus("For rrblup method, the example results as follows:(Better to add some discription about the example data.)")
      
      output$exampleplot1 <- renderImage({
        list(src = paste0("./example/rrblup/result/traits.png"), contentType = "image/png", height = 400)
      }, deleteFile = FALSE)
      
      output$exampleplot2 <- renderImage({
        req(input$rrbluptrait)
        list(src = paste0("./example/rrblup/result/plot_", input$rrbluptrait, ".png"), contentType = "image/png", height = 400)
      }, deleteFile = FALSE)
      
      output$exampleplot3 <- renderImage({
        req(input$rrbluptrait)
        list(src = paste0("./example/rrblup/result/plot2_", input$rrbluptrait, ".png"), contentType = "image/png", height = 400)
      }, deleteFile = FALSE)
      
      EDataloc("The discription of the plots.")
      
    } else if (input$model == "DNNGP3") {
      runStatus("For DNNGP method, the example results as follows:(Better to add some discription about the example data.)")
      
      output$exampleplot1 <- renderImage({
        list(src = paste0("./example/rrblup/result/traits.png"), contentType = "image/png", height = 400)
      }, deleteFile = FALSE)
    } else if (input$model == "DeepGS") {
      
    } else {
      
    }
  })
  
  session$onSessionEnded(function() {
    stop("Session ended")
  }) 
  
}
shinyApp(ui = ui, server = server)
