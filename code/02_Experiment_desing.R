# libraries ####
source("code/00_packages.R")
source("code/global_variables.R")

# Experiment design ####
# You just need modify script in Experiment heading

label <-    c("lfq_wt_1", "lfq_wt_2", "lfq_wt_3", "lfq_wt_4",
              "lfq_cln3_lux1_1", "lfq_cln3_lux1_2", "lfq_cln3_lux1_4",
              "lfq_cln3_lux2_2", "lfq_cln3_lux2_3", "lfq_cln3_lux2_4",
              "lfq_cln12_2", "lfq_cln12_3",  "lfq_cln12_4")

label_columns <-  label

columns_to_rename <- c("WT_1","WT_2","WT_3", "WT_4",
                       "CLN3_Lux1_1", "CLN3_Lux1_2", "CLN3_Lux1_4",
                       "CLN3_Lux2_2", "CLN3_Lux2_3", "CLN3_Lux2_4",
                       "CLN12_2", "CLN12_3", "CLN12_4")


condition <- as.factor(c("WT","WT","WT","WT",
                         "CLN3_Lux1","CLN3_Lux1","CLN3_Lux1",
                         "CLN3_Lux2","CLN3_Lux2","CLN3_Lux2",
                         "CLN12","CLN12","CLN12"))
replicate <- c(1,2,3,4,
               1,2,4,
               2,3,4,
               2,3,4)

experiment <- c("Brains","Brains","Brains","Brains",
                "Brains","Brains","Brains",
                "Brains","Brains","Brains",
                "Brains","Brains","Brains")

Exp_design <- cbind.data.frame(label,condition,replicate, experiment)
Exp_design

print(Exp_design)

Exp_design %>% 
  write_xlsx(paste0(output_path,"tables/experiment_desing.xlsx"))
