# Purpose: Infer potential ligand-receptor communication with optional CellChat.
# Input:   results/objects/04_malignancy_candidates.rds
# Output:  CellChat object, interaction table, and network plot
# Caution: Results are computational predictions, not direct interaction evidence.

source("functions/workflow_utils.R")
source("functions/plotting_functions.R")
cfg <- read_config()
if (!isTRUE(cfg$communication$enabled)) stop("Optional module disabled. Set communication.enabled: true after installing CellChat.", call. = FALSE)
assert_packages(c("CellChat", "Seurat"), optional = TRUE)
object <- load_seurat("04_malignancy_candidates", cfg)
validate_metadata(object, cfg$communication$group_by)
data_input <- Seurat::GetAssayData(object, assay = "RNA", layer = "data")
meta <- object[[]]
cellchat <- CellChat::createCellChat(data_input, meta = meta, group.by = cfg$communication$group_by)
database <- if (tolower(cfg$project$species) == "human") CellChat::CellChatDB.human else CellChat::CellChatDB.mouse
cellchat@DB <- CellChat::subsetDB(database, search = cfg$communication$database)
cellchat <- CellChat::subsetData(cellchat)
cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
cellchat <- CellChat::computeCommunProb(cellchat)
cellchat <- CellChat::filterCommunication(cellchat, min.cells = cfg$communication$min_cells)
cellchat <- CellChat::computeCommunProbPathway(cellchat)
cellchat <- CellChat::aggregateNet(cellchat)
saveRDS(cellchat, project_path(cfg$project$output_dir, "objects", "07_cellchat.rds"))
write_table(CellChat::subsetCommunication(cellchat), "07_predicted_communications.csv", cfg)
grDevices::pdf(project_path(cfg$project$output_dir, "figures", "07_cellchat_network.pdf"), width = 6, height = 6)
CellChat::netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)),
                          weight.scale = TRUE, label.edge = FALSE)
grDevices::dev.off()
save_session_info("07_Cell_Cell_Communication", cfg)
