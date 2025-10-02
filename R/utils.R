exctract_intersect_metabolites <- function(cleaned_labels, reference) {
  # Use Levenshtein similarity to get similarities between pre-processed strings
  # (removed special characters)
  # Afterwards a manual comparison between non-100% similar metabolites
  
  # load uk biobank metabolites names from supplementary material
  SuppInfo_Table1_metabolic_variables_description <- read_excel(reference)
  
  
  # only execute if the matchframe csv does not exist
  if (!file.exists("Matchframe_backup.csv")) {
    library(RecordLinkage) # to compare strings
    # compare name of metabolites with those in supp. material
    # levenshtein distance
    compare_str_vector <- function(Vstr1, Vstr2) {
      n = length(Vstr1) # metabolites
      p = length(Vstr2) # suppl names
      similarities = matrix(0, nrow = n, ncol = p)
      for (i in 1:n) {
        for (j in 1:p) {
          word1 = Vstr1[i] %>% tolower() %>% str_remove_all("[:.-]")
          word2 = Vstr2[j] %>% tolower()  %>% str_remove_all("[:.-]")
          
          # Remove (mol/l) from both strings
          word1 = word1 %>% str_remove_all("\\(mol/l\\)") %>% str_remove_all("\\(mmol/l\\)") %>% str_remove_all("mmol/l") 
          word2 = word2 %>% str_remove_all("\\(mol/l\\)") %>% str_remove_all("\\(mmol/l\\)") %>% str_remove_all("mmol/l") 
          
          # since "ratio" is used in word1, where in word2 it is percentage, we change
          # Handle ratio and percentage in word1
          if (grepl("ratio", word1)) {
            # Remove "ratio" and add "%" at the end
            word1 = str_remove(word1, "ratio") %>% paste0("%")
          } else if (grepl("percentage", word1)) {
            # Remove "percentage" and add "%" at the end
            word1 = str_remove(word1, "percentage") %>% paste0("%")
          }
          
          # Handle ratio and percentage in word2
          if (grepl("ratio", word2)) {
            # Remove "ratio" and add "%" at the end
            word2 = str_remove(word2, "ratio") %>% paste0(" %")
          } else if (grepl("percentage", word2)) {
            # Remove "percentage" and add "%" at the end
            word2 = str_remove(word2, "percentage") %>% paste0(" %")
          }
          
          # sort the words of both strings alphabetically
          word1 = word1 %>% scan(text = ., what = " ", quiet = T)# %>% sort()
          word2 = word2 %>% scan(text = ., what = " ", quiet = T)# %>% sort()
          

          if ("total" %in% word1 & "free" %in% word2 || "free" %in% word1 & "total" %in% word2) {
            similarities[i,j] = 0 # free and total are different measurements
          } else {
            word1 = word1 %>% paste(collapse = " ") %>% str_remove_all("[ ]")
            word2 = word2 %>% paste(collapse = " ") %>% str_remove_all("[ ]")
            similarities[i,j] = levenshteinSim(word1, word2)
          }
          
          #wordsimmat = matrix(NA, nrow = length(word1), ncol = length(word2))
          ## get mean levenschtein sim for all words
          #for (wi in 1:length(word1)) {
          #  for (wj in 1:length(word2)) {
          #    w1 = testword[1:wi] %>% paste(collapse = " ")
          #    w2 = testword2[1:wj] %>% paste(collapse = " ")
          #    wordsimmat[wi,wj] = levenshteinSim(w1,w2)
          #  }
          #}
          
        }
      }
      colnames(similarities) = Vstr2
      rownames(similarities) = Vstr1
      return(similarities)
    }
    
    # Get matrix of string similarities
    str_sims = compare_str_vector(cleaned_labels, SuppInfo_Table1_metabolic_variables_description$Description)
    
    # Dataframe with matches
    n_data_metabolites = length(cleaned_labels)
    matchframe = data.frame(data_metabolite = rep(NA, n_data_metabolites), validation_metabolite = rep(NA, n_data_metabolites), similarity = rep(0, n_data_metabolites))
    
    # for each row, take the highest match (may not always be identical)
    for (m in 1:n_data_metabolites) {
      name = cleaned_labels[m]
      match = which.max(str_sims[m,])
      value = str_sims[m,match]
      match = match %>% names # extract column name
      matchframe[m,] = c(name, match, value)
    }
    matchframe$similarity = as.double(matchframe$similarity)
    
    # Manually assign whether the metabolites are true matches
    matchframe["true.match"] = NA
    for (i in 1:nrow(matchframe)) {
      # No need for comparison of 100% similarity
      if (matchframe[i,3] == 1) {
        matchframe[i, "true.match"] = 1
      } else if (matchframe[i,3] == 0) { # Skip when 0% match
        matchframe[i, "true.match"] = 0
      } else {
        # Display the current row for manual check
        cat("Row", i, "values:\n")
        cat("Column 1:", as.character(matchframe[i,1]), "\n") 
        cat("Column 2:", as.character(matchframe[i,2]), "\n")
        cat("Similarity value:", as.character(matchframe[i,3]), "\n")
        
        user_response = readline(prompt = "Correct (False: 0 / True: 1 / Doubt: 2)? ")
        
        matchframe[i, "true.match"] = user_response
        
        cat("\014")  # This clears the console for better overview
      }
    }
    # remove the "NULL"row
    matchframe = matchframe[1:229,]
    
    # Rename column 2 from the suppinfo to match the column name in matchframe
    colnames(SuppInfo_Table1_metabolic_variables_description)[2] = "validation_metabolite"
    
    # add metabolic var to new column from suppinfo
    matchframe <- matchframe %>%
      left_join(
        SuppInfo_Table1_metabolic_variables_description %>% select(validation_metabolite, `Metabolic variable`), 
        by = "validation_metabolite"
      )
    # Add a column with the metabolic names in the data for the descriptions
    matchframe$column_name = pre_meal_metab %>% names()
    matchframe$column_name = gsub("_.*$", "", matchframe$column_name)
    
    
    # save matchfame to csv (backup)
    write.csv(matchframe, "Matchframe_backup.csv", )
  } else { # load to avoid re-processing
    matchframe = read.csv("Matchframe_backup.csv", header = T)
  }
  
  # we now have a mask for match frame
  matchframe$true.match = as.integer(matchframe$true.match)
  matchframe %>% filter(true.match == 1) %>% head
  
  # number of identical metabolites between the two
  matchframe %>% filter(true.match == 1) %>% nrow
  
  # total metabolites in ukbiobank model (survRRR)
  SuppInfo_Table1_metabolic_variables_description %>% nrow
  
  # ratio of total in biobank model
  (matchframe %>% filter(true.match == 1) %>% nrow) / 
    (SuppInfo_Table1_metabolic_variables_description %>% nrow)
  return(matchframe)
}


find_optimal_exclusion <- function(dataset) {
  # Thresholds at which, when a person misses more than thresh bio marker measurements,
  # the patient is excluded
  # UNUSED FUNCTION. USED IMPUTATION METHOD MIN(c)/2
  exclustion_thresh = seq(0, 1, 0.001)
  
  # Dimensions of the original set -> which to optimize on
  og_dim = dim(dataset)
  p = og_dim[2] - 1
  p_rows = rep(NA, length(exclustion_thresh))
  p_cols = p_rows
  i = 1
  
  for (thresh in exclustion_thresh) {
    # Remove patients with more than thesh missing measurements for bio markers
    dataset_i <- dataset[!(((rowSums(dataset %>% is.na()))/p) > thresh),]
    
    # Remove metabolites with more than 1% missing measurements
    dataset_i <- dataset_i[,!(((colSums(dataset_i %>% is.na()))/nrow(dataset)) > 0.01)]
    
    # Remove people with missing metabolite values
    dataset_i <- dataset_i[!(((rowSums(dataset_i %>% is.na()))/p) > 0.00),]
    
    p_rows[i] = nrow(dataset_i) / og_dim[1]
    p_cols[i] = (ncol(dataset_i) - 1) / p
    i = i + 1
  }
  
  res_df = data.frame(`Exclusion.threshold` = exclustion_thresh, Rows = p_rows, Columns = p_cols)
  
  p = pivot_longer(res_df, cols = c(Rows, Columns)) %>%
    mutate(`Proportion kept` = value) %>%
    mutate(`dimension` = name) %>%
    ggplot(aes(x = `Exclusion.threshold`, y = `Proportion kept`, col = `dimension`)) +
    geom_line() +
    theme_minimal()
  print(p)
  return(res_df)
}


make_long <- function(occasion = 1, data, outcome_death, outcome_MI,
                      outcome_TIA, meal_covs, standardize = FALSE, filename = "dlong", log_transform = F) {
  ## occasion in {0, 1} -- 1: Pre-meal, 2: Post-meal
  ## data               -- full data.frame
  ## outcome_death      -- extracted death times
  ## outcome_MI         -- Extracted Myocardial infarction info
  ## outcome_TIA        -- Extracted Transient ischemic attack info
  ## meal_covs          -- Metabolite covariates related to prandial state
  ## standardize {T/F}  -- Whether to standardize the metabolites
  ## filename           -- Name of the csv file to write it to {name}.csv
  ## log_transform      -- Whether to log transform before scale and center
  
  # get the file prefix related to occasion
  prefix = c("pre", "post")[occasion]
  # Add to filename whether the covariates are standardized
  sdized = c("", "_standardized")[standardize %>% as.integer() + 1]
  
  # if dlong does not exist, create it 
  if (!file.exists(glue("//vf-DataSafe/DataSafe$/Div0/ict/survRRR in NEO_2873/{filename}.csv"))) {
    # correctly scale the covariates
    cov_descr = meal_covs %>% lapply(function(x) attributes(x)$label) %>% as.character()
    for (m_i in 1:(length(cov_descr)-1)) { # -1 due to id column at the end
      # avoid changing ratio columns 
      if (TRUE) { #(!grepl("%", cov_descr[m_i], fixed = T)) {
        process_column = meal_covs[m_i]
        if (log_transform) {
          # add 1
          # log transform
          process_column = process_column + 1
          process_column = log(process_column) # no log for survRRR
        }
        
        # center and scale
        process_column = scale(process_column, center = T, scale = T)
        meal_covs[,m_i] = process_column # De_COMMENT THIS!
      }
    }
    
    # Create empty dataframe to store results
    dlong = data.frame(id = numeric(), Tstart = numeric(), Tstop = numeric(),
                       status = numeric(), from = numeric(), to = numeric(),
                       trans = numeric())
    
    # Add metabolite columns as covariates to dlong
    metab_columns <- setdiff(names(meal_covs), "id")
    for (col in metab_columns) {
      dlong[[col]] <- numeric()
    }
    
    # For each individual in meal_covs we assemble the long dataframe
    for (id_i in meal_covs$id) {
      # print(glue("Processing id {id_i} out of {length(meal_covs$id)}"))
      # Get outcome data for current ID
      visit_id = data %>% filter(id == id_i)
      
      # Skip if no matching record found
      # if (nrow(visit_id) == 0) {
      #   next
      # }
      
      # Skip record of no final date
      if (is.na(visit_id$einddatum2)) {
        next
      }
      
      # Age at visit
      Tstart = visit_id$leeftijd
      
      # Get patient's metabolite values
      patient_metabs <- meal_covs %>% filter(id == id_i)
      
      for (trans in 1:5){ #outcomes: DIA, DTH, MI, TIA, FIRST OUTCOME
        # Set default from state (1 = healthy)
        from_state = 1
        # Set destination state based on transition
        to_state = trans + 1  # DIA=2, DTH=3, MI=4, TIA=5, first = 6
        
        if (trans == 1) { #DIA
          if (visit_id$diab_prev != 1) { # Only get outcome if patient has no diabetes before study
            status = visit_id$diabetes2
            Tstop =  as.numeric(# Years passed since start of study
              difftime(visit_id$diabetes2_date, visit_id$visitdd, units = "days")
            )/365.25 # For leap years
            if (is.na(Tstop)) { # Censored
              # Years passed since start of study
              Tstop =  as.numeric(
                difftime(visit_id$einddatum2,visit_id$visitdd,units = "days")
              )/365.25
            } 
          } # Else avoid diabetes for the patient (before study) 
          else {
            next
          }
        }
        
        if (trans == 2) { # DeaTH
          # Retrieve info whether person has died
          status = outcome_death %>% filter(id == id_i) %>%
            select(outcome_DTH) %>% as.numeric() # 0 = not censored = event happened
          if (status == 1) {
            # Extract death date
            Tstop = outcome_death %>% filter(id == id_i) %>%
              pull(outcome_DTH_date) %>% as.Date()
            Tstop =  as.numeric(# Years passed since start of study
              difftime(Tstop,visit_id$visitdd,  units = "days")
            )/365.25
          } else { # Censored
            # Years passed since start of study
            Tstop =  as.numeric(
              difftime(visit_id$einddatum2,visit_id$visitdd,units = "days")
            )/365.25
          }
        }
        
        if (trans == 3) { # MI
          #status = outcome_MI %>% filter(id == id_i) %>%
          #  select(MI2_inc) %>% as.numeric()
          
          # If date1 is not empty (not NA), the person had MI
          status = ifelse(!is.na(outcome_MI %>% filter(id == id_i) %>% pull(MI2_Date_1)), 1, 0)
          if (status == 1) {
            Tstop = outcome_MI %>% filter(id == id_i) %>%
              pull(MI2_Date_1) %>% as.Date()
            Tstop =  as.numeric(# Years passed since start of study
              difftime(Tstop,visit_id$visitdd,  units = "days")
            )/365.25
            
            # If the person had MI before study, skip this outcome for this person
            if (Tstop < 0) {
              next
            }
          } else { # Censored
            # Years passed since start of study
            Tstop =  as.numeric(
              difftime(visit_id$einddatum2,visit_id$visitdd,units = "days")
            )/365.25
          }
        }
        
        if (trans == 4) { # TIA
          # status = outcome_TIA %>% filter(id == id_i) %>%
          #   select(TIA2_inc) %>% as.numeric()
          
          # If date1 is not empty, the person had TIA
          status = ifelse(!is.na(outcome_TIA %>% filter(id == id_i) %>% pull(TIA2_date_1)), 1, 0)
          
          
          if (status == 1) {
            Tstop = outcome_TIA %>% filter(id == id_i) %>%
              pull(TIA2_date_1) %>% as.Date()
            Tstop =  as.numeric(# Years passed since start of study
              difftime(Tstop,visit_id$visitdd,  units = "days")
            )/365.25
            # If the date < entry date, skip (person had it before study)
            if (Tstop < 0) {
              next
            }
          } else { # Censored
            # Years passed since start of study
            Tstop =  as.numeric(
              difftime(visit_id$einddatum2,visit_id$visitdd,units = "days")
            )/365.25
          }
          
        }
        
        if (trans == 5) { # First
          # Get all events for patient id_i from dlong and choose whichever happens first
          new_row = dlong %>% filter(id == id_i) # Events of id_i
          
          # Pass if person has no measurements
          if (nrow(new_row) == 0) {
            next
          }
          
          # Check if patient had an event [something goes wrong here, I think]
          if (sum(new_row$status) > 0) { 
            # Pick the event that occurred earliest
            new_row = new_row %>% filter(status == 1) %>% slice_min(Tstop, n = 1)
            # Ensure that it is just one[1] row
            new_row = new_row[1,]
          } else {
            # Censored, just take any row
            new_row = new_row %>% slice_max(Tstop, n = 1, na_rm = TRUE)
            # Ensure that it is just one[1] row
            new_row = new_row[1,]
          }
          
          
          # Store specifics from this transition into that first event
          new_row$from = from_state
          new_row$to = to_state
          new_row$trans = trans
          
          # We just steal the metabolic variables from the first occurrence
          
        } else { # Have to create a new row from scratch
          # Create new row
          new_row <- data.frame(
            id = id_i,
            Tstart = Tstart,  # Age of person when included
            Tstop = Tstart + Tstop, # Age of person at stop
            status = status,
            from = from_state,
            to = to_state,
            trans = trans
          )
          # Add metabolite values for this patient
          if (nrow(patient_metabs) > 0) {
            for (col in metab_columns) {
              new_row[[col]] <- patient_metabs[[col]][1]
            }
          }
        }
        
        
        
        
        
        # Add patient metabolic covariates
        # Add the row to dlong
        dlong <- rbind(dlong, new_row)
        
        # if (id_i >= 500) {
        #   test_dlong <<- dlong
        #   stop("Test your dlong")
        # }
        
      }
    }
    # Write it to csv in secure location
    write.csv(dlong, glue("//vf-DataSafe/DataSafe$/Div0/ict/survRRR in NEO_2873/{filename}.csv"))
  } else { # Not needed, too much trouble to load (why did I write this?)
    dlong = read.csv(glue("//vf-DataSafe/DataSafe$/Div0/ict/survRRR in NEO_2873/{filename}.csv"), header = T)[,-1]
  }
  return(dlong)
}


calculate_cindex <- function(ext_data, transition, coef_table, outcome_map) {
  # Function to calculate C-index for a specific transition
  # ext_data:     data for external validation (dataframe)
  # transition:   index of transition (outcome)
  # coef_table:   alphas from ukbiobank
  # outcome_map:  string values of transition indexes
  
  # Filter data for the specific transition (separate dataset)
  transition_data = ext_data %>% filter(trans == transition)
  
  # Get outcome name from the mapping
  outcome_name = outcome_map[transition]
  
  # Get coefficients for this outcome and rename them according to external dataset
  outcome_coefs = coef_table$coefficient
  names(outcome_coefs) = coef_table$metabolic_variable
  
  # Grab the metabolites names from the external dataset
  metabolic_vars = setdiff(names(transition_data), 
                           c("status", "Tstart", "Tstop", "trans", 
                             "from", "to", "id", "X"))
  
  
  # Filter coefficients to only use available metabolic variables, also reorder in order of external
  usable_coefs = outcome_coefs[metabolic_vars]
  print(length(usable_coefs))
  usable_coefs = na.omit(usable_coefs)
  print(length(usable_coefs))
  
  # Matrix multiplication to get linear predictor, then convert to vector
  linear_pred = as.vector(as.matrix(transition_data[, names(usable_coefs)]) %*% usable_coefs)
  print(length(linear_pred))
  
  # Calculate the linear predictors
  transition_data$p_lps = linear_pred
  transition_data$n_lps = -linear_pred
  
  
  ## Uno's C
  # Survival object with follow-up time due to limitations in Uno's C
  Uno_surv_obj = Surv(transition_data$Tstop - transition_data$Tstart, transition_data$status)
  
  # Compute Uno's C (n/G2)
  cpos_uno = concordance(Uno_surv_obj ~ p_lps, timewt="n/G2", data = transition_data)
  cneg_uno = concordance(Uno_surv_obj ~ n_lps, timewt="n/G2", data = transition_data)
  
  ## Harrells C
  # Survival object with age scale
  Har_surv_obj = Surv(transition_data$Tstart, transition_data$Tstop, transition_data$status)
  
  cpos_har = concordance(Har_surv_obj ~ p_lps, data = transition_data)
  cneg_har = concordance(Har_surv_obj ~ n_lps, data = transition_data)
  
  ## Harrels C follow-up time
  # Survival object with follow-up time for comparison
  Har_ft_surv_obj = Surv(transition_data$Tstop - transition_data$Tstart, transition_data$status)
  
  cpos_har_ft = concordance(Har_ft_surv_obj ~ p_lps, data = transition_data)
  cneg_har_ft = concordance(Har_ft_surv_obj ~ n_lps, data = transition_data)
  
  Cframe = list()
  
  # Determine which direction is better
  # for (Cm in c("Uno (Follow-up)", "Harrell (Age)", "Harrell (Follow-up)")) {
  print("Uno (Follow-up)")
  if (cpos_uno$concordance >= cneg_uno$concordance) {
    cat(glue("For outcome '{outcome_name}', used positive direction (C-index = {round(cpos_uno$concordance, 4)})."), "\n")
    Cframe = append(Cframe, cpos_uno)
    print(cpos_uno)
  } else {
    cat(glue("For outcome '{outcome_name}', used negative direction because C-index < 0.5 (original = {round(cpos_uno$concordance, 4)}, flipped = {round(cneg_uno$concordance, 4)})."),"\n")
    Cframe = append(Cframe, cneg_uno)
    print(cneg_uno)
  }
  
  print("Harrell (Age)")
  if (cpos_har$concordance >= cneg_har$concordance) {
    cat(glue("For outcome '{outcome_name}', used positive direction (C-index = {round(cpos_har$concordance, 4)})."), "\n")
    Cframe = append(Cframe, cpos_har)
    print(cpos_har)
  } else {
    cat(glue("For outcome '{outcome_name}', used negative direction because C-index < 0.5 (original = {round(cpos_har$concordance, 4)}, flipped = {round(cneg_har$concordance, 4)})."),"\n")
    Cframe = append(Cframe, cneg_har)
    print(cneg_har)
  }
  
  print("Harrel (Follow-up)")
  if (cpos_har_ft$concordance >= cneg_har_ft$concordance) {
    cat(glue("For outcome '{outcome_name}', used positive direction (C-index = {round(cpos_har_ft$concordance, 4)})."), "\n")
    Cframe = append(Cframe, cpos_har_ft)
    print(cpos_har_ft)
  } else {
    cat(glue("For outcome '{outcome_name}', used negative direction because C-index < 0.5 (original = {round(cpos_har_ft$concordance, 4)}, flipped = {round(cneg_har_ft$concordance, 4)})."),"\n")
    Cframe = append(Cframe, cneg_har_ft)
    print(cneg_har_ft)
  }
  return(Cframe)
}