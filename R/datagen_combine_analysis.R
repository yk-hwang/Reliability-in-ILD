


library(MplusAutomation)
library(dplyr);  library(stringr); library(purrr); library(readr); library(openxlsx)


base_dir   <- "C:/Users/hwang/OneDrive/Desktop/sim8"

sets<- 1:5

logv1_vals <- c(4.76190,1.78571,0.97402,0.06252,0.13227)
logv2_vals <- c(1,1,1,1,1)


template_code <- "
TITLE: mear(1) new condition

MONTECARLO:
  NAMES = y; 
  NREPS = 100; 
  NOBS = 300; 
  NCSIZES = 1; 
  SEED = 12345; 
  CSIZES = 20(15); 
  REPSAVE = ALL; 
  SAVE = mearsim_cond12.rep*.dat; 

ANALYSIS: TYPE = TWOLEVEL RANDOM; ESTIMATOR = BAYES; PROCESSORS = 2; FBITER = (120000);

MODEL POPULATION:
  %%WITHIN%%
    f BY y@1(&1);
    phi | f ON f&1;
    [f@0];
    
    y*%s; !measurement error
    f*%s; !innovation 
    
  %%BETWEEN%%
    [y*6.5]; 
    y*3;

    [phi*.4];
    phi*0.000001;          ! random effect

    ! no correlation
    y phi WITH y*0 phi*0;
    
    MODEL: 
 %%WITHIN%%
    f BY y@1(&1);
    phi | f ON f&1;
    [f@0];
    
    y; !measurement error
    f; !innovation 
    
  %%BETWEEN%%
    [y]; 
    y;

    [phi];
    phi;          ! random effect

    ! no correlation
    y phi WITH y*0 phi*0;
    
"

old_wd <- getwd()
for (i in seq_along(sets)) {
  set_i <- sets[i]
  dir_i <- file.path(base_dir, paste0("set", set_i))
  inp_i <- file.path(dir_i, "mearsim_cond12.inp")
  
  if (!dir.exists(dir_i)) dir.create(dir_i, recursive = TRUE)
  
  mplus_code <- sprintf(
    template_code,
    logv1_vals[i],logv2_vals[i]
  )
  

  writeLines(mplus_code, inp_i)
  setwd(dir_i)
  runModels(target = basename(inp_i))
  setwd(old_wd)
}
message("▶ 5세트 Monte Carlo 실행 완료")

dirs         <- file.path(base_dir, paste0("set", sets))
combined_dir <- file.path(base_dir, "combined")
if (!dir.exists(combined_dir)) dir.create(combined_dir)

combined_data <- vector("list", length = 100)
names(combined_data) <- paste0("rep", 1:100)

for (rep in seq_along(combined_data)) {
  dfs <- setNames(
    lapply(dirs, function(d) {
      fn <- file.path(d, paste0("mearsim_cond12.rep", rep, ".dat"))
      read.table(fn, header = FALSE)
    }),
    paste0("set", sets)
  )
  
  stopifnot(names(dfs) == paste0("set", sets))
  
  
  for (set_nm in names(dfs)) {
    offset <- (as.numeric(sub("set", "", set_nm)) - 1) * 20
    dfs[[set_nm]][, 2] <- dfs[[set_nm]][, 2] + offset
  }
  

  clusters <- sort(unique(unlist(lapply(dfs, `[[`, 2))))
  

  ordered_list <- vector("list", length(clusters) * length(dfs))
  idx <- 1
  for (cl in clusters) {
    for (set_nm in names(dfs)) {
      block <- dfs[[set_nm]][ dfs[[set_nm]][,2] == cl, , drop = FALSE ]
      if (nrow(block) > 0) {
        ordered_list[[idx]] <- block
        idx <- idx + 1
      }
    }
  }
  ordered_list <- ordered_list[seq_len(idx-1)]
  

  combined_df <- do.call(rbind, ordered_list)
  combined_data[[rep]] <- combined_df
  write.table(
    combined_df,
    file      = file.path(combined_dir, paste0("mearsim_cond12.rep", rep, ".dat")),
    row.names = FALSE,
    col.names = FALSE,
    quote     = FALSE
  )
}

message("▶ 5세트 결과 병합 완료: cluster×set 순서로 정렬된 crp*.dat 생성")



setwd("C:/Users/hwang/OneDrive/Desktop/sim8/combined/")
files <- list.files(pattern="^mearsim_cond12\\.rep\\d+\\.dat$"); 
files <- files[order(as.integer(sub(".*rep(\\d+)\\.dat$","\\1",files)))]
writeLines(files, "mearsim_cond12.replist.dat")

library(MplusAutomation)
library(readr)
library(stringr)
library(dplyr)
library(purrr)
library(openxlsx)

mplus_monte_carlo_code <- "
TITLE: mear(1) condition9

data: 
file is mearsim_cond12.replist.dat; 
type = montecarlo; 

variable: NAMES = y id; cluster = id; 

ANALYSIS: TYPE = TWOLEVEL RANDOM; ESTIMATOR = BAYES; PROCESSORS = 2; FBITER = (120000);

MODEL:
  %WITHIN%
  f BY y@1(&1);
  phi | f ON f&1;
  [f@0];

  logv1 | y; 
  logv2 | f;

  %BETWEEN%
  [y]; 
   y;             
  [phi*.4];
   phi;             ! random effect
  [logv1];
   logv1;           ! measurement err var cluster mean
  [logv2];
   logv2;           ! innovation var cluster mean

  ! no correlation
  y phi logv1 logv2 WITH y*0 phi*0 logv1*0 logv2*0;
"


writeLines(mplus_monte_carlo_code, con = "mearsim_cond12.inp")
runModels(target = "mearsim_cond12.inp")
var_names <- c(
  "Y                                F10.3",
  "F Mean                           F10.3",
  "F Median                         F10.3",
  "F Standard Deviation             F10.3",
  "F 2.5% Value                     F10.3",
  "F 97.5% Value                    F10.3",
  "F&1 Mean                         F10.3",
  "F&1 Median                       F10.3",
  "F&1 Standard Deviation           F10.3",
  "F&1 2.5% Value                   F10.3",
  "F&1 97.5% Value                  F10.3",
  "PHI Mean                         F10.3",
  "PHI Median                       F10.3",
  "PHI Standard Deviation           F10.3",
  "PHI 2.5% Value                   F10.3",
  "PHI 97.5% Value                  F10.3",
  "LOGV1 Mean                       F10.3",
  "LOGV1 Median                     F10.3",
  "LOGV1 Standard Deviation         F10.3",
  "LOGV1 2.5% Value                 F10.3",
  "LOGV1 97.5% Value                F10.3",
  "LOGV2 Mean                       F10.3",
  "LOGV2 Median                     F10.3",
  "LOGV2 Standard Deviation         F10.3",
  "LOGV2 2.5% Value                 F10.3",
  "LOGV2 97.5% Value                F10.3",
  "B_Y Mean                         F10.3",
  "B_Y Median                       F10.3",
  "B_Y Standard Deviation           F10.3",
  "B_Y 2.5% Value                   F10.3",
  "B_Y 97.5% Value                  F10.3",
  "ID                               I4",
  "_TRUETIME                        I3",
  "_TIMEPOINT                       I3"
)


var_names <- var_names %>%
  str_replace_all("F10.3|I4|I3", "") %>%  
  str_trim() %>%                        
  str_replace_all("\\s+", " ") %>%       
  str_replace_all(" ", "")              

reps <- 1:100
varnames_required <- c('PHIMean', 'LOGV1Mean', 'LOGV2Mean', 'ID')
process_data <- function(rep, estimates) {
  data_file <- paste0("test_", rep, ".dat")
  if (!file.exists(data_file)) return(NULL)
  
  raw_df <- read.table(data_file, header = FALSE) 
  
  names(raw_df) <- var_names                                                                                                                                                                                                                                                                                                                                                                     
  
  df_filtered <- raw_df %>% select(all_of(varnames_required))
  raw_n <- nrow(df_filtered)
  df <- df_filtered %>% filter(PHIMean < 1)
  removed_count <- raw_n - nrow(df)

  df <- df %>%
    mutate(
      within_variance = exp(estimates$meanslogv2 + LOGV2Mean) / (1 - PHIMean^2),
      innovv = exp(estimates$meanslogv2 + LOGV2Mean),
      meav = exp(estimates$meanslogv1 + LOGV1Mean),
      rel_w = within_variance / (within_variance + meav)
    ) %>% mutate(rel_w_variance = var(rel_w)) 

  df_m <- df %>%
    select(-ID) %>%
    summarise(across(everything(), mean)) %>%
    mutate(
      relb = estimates$variancey / (estimates$variancey + within_variance + meav),
      rep = rep,
      removed_subjects = removed_count
    ) %>% mutate(relb_variance = var(relb))
  
  return(list(summary = df_m, details = df))
}

main_func <- function(rep) {
  input_file <- sprintf("mearsim_cond12.rep%d.inp", rep)
  data_file  <- sprintf("mearsim_cond12.rep%d.dat", rep)
  output_file <- sprintf("mearsim_cond12.rep%d.out", rep)
  
  if (!file.exists(data_file)) return(NULL)
  
  mplus_code <- paste0(
    "TITLE: mear(1) analysis\n",
    "DATA: FILE = ", data_file, ";\n",
    "VARIABLE: NAMES = y id; CLUSTER = id;\n",
    "ANALYSIS: TYPE = TWOLEVEL RANDOM; ESTIMATOR = BAYES; PROCESSORS = 2; FBITER = (120000);\n",
    "MODEL:\n",
    "%WITHIN%\n",
    "f BY y@1 (&1);\n",
    "phi | f ON f&1;\n",
    "[f@0];\n",
    "logv1 | y;\n",
    "logv2 | f;\n",
    "%BETWEEN%\n",
    "[y]; y; [phi]; phi; [logv1]; logv1; [logv2]; logv2;\n",
    "y phi logv1 logv2 WITH y*0 phi*0 logv1*0 logv2*0;\n",
    "OUTPUT: RESIDUAL(CLUSTER);\n",
    
    "SAVEDATA: FILE IS test_", rep, ".dat; SAVE = fscores(1);\n"
  )
  
  writeLines(mplus_code, con = input_file)
  runModels(target = input_file)
  
  if (!file.exists(output_file)) return(NULL)
  
  out <- readModels(output_file)$parameters$unstandardized
  meanslogv1 <- out$est[out$paramHeader=="Means" & out$param=="LOGV1"]
  meanslogv2 <- out$est[out$paramHeader=="Means" & out$param=="LOGV2"]
  variancey  <- out$est[out$paramHeader=="Variances" & out$param=="Y"]
  
  estimates <- list(meanslogv1=meanslogv1, meanslogv2=meanslogv2, variancey=variancey)
  
  process_data(rep, estimates)
}


res_list <- map(reps, main_func)

final_results <- map_df(res_list, "summary")


#write.xlsx(final_results, "reliability.xlsx", overwrite = TRUE)

#if (file.exists("reliability.xlsx")) {
#  old_results <- read.xlsx("reliability.xlsx")
#  final_results <- bind_rows(old_results, final_results)
#}
#write.xlsx(final_results, "reliability.xlsx")
mean_reliability <- mutate(final_results, mean_relb = mean(relb), var_relb = mean(relb_variance), 
                           mean_relw = mean(rel_w), var_relw = mean(rel_w_variance),
                           over15=length(which(removed_subjects > 15)),
                           same15=length(which(removed_subjects == 15)),
                           zero=length(which(removed_subjects == 0)),
                           mean_phi = mean(PHIMean),
                           var_phi = var(PHIMean))
write.xlsx(mean_reliability, "mean_reliability.xlsx", overwrite = TRUE)
#View(res_list[[3]]$details)


