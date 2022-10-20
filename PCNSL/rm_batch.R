Sys.setenv(LANG = "en") # change language
library(sva)

## (1) read data in prospective cohort and retrospective cohort
pc <- read.csv("./data/batch/PC_normalized.csv")
rc <- read.csv("./data/batch/RC_normalized.csv")

###PFS and label
pc_pfs <- read.csv("./data/batch/prospect_full.csv", row.names = 1)
pc_pfs <- pc_pfs[1:2, ]
rc_pfs <- read.csv("./data/batch/retro_full.csv", row.names = 1)
rc_pfs <- rc_pfs[1:2, ]


#merge dataframe
total <- merge(x = pc, y = rc, by = "metabolite", all = FALSE)
rownames(total) <- total$metabolite

if ("metabolite" %in% colnames(total)){
    total <- total[, -1]
}else{
    print("!!")
}

### (2)remove batch effect
## read groups
group_df <- read.csv("./data/batch/group.csv")
batch <- group_df$Batch

## other covariate
tissue <- group_df$Group
table(batch,tissue)
modCombat <- model.matrix(~tissue)

##remove batch effect
expr_batch <- ComBat(dat = total, batch = batch, mod = modCombat)

rc_sample <- group_df[group_df$Batch == 1,  "Sample"]
pc_sample <- group_df[group_df$Batch == 2,  "Sample"]


rc_batch <- rbind(rc_pfs , expr_batch[, rc_sample])
pc_batch <- rbind(pc_pfs, expr_batch[, pc_sample])

write.csv(x = rc_batch, file = "./data/batch/rc_rm_batch.csv")
write.csv(x = pc_batch, file = "./data/batch/pc_rm_batch.csv")
