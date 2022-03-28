# ----------------------------------------------------------------------------
# create_tree_labels.R 
# Author: Antoine Beauchamp
# Created: February 16th, 2021
#
# Description
# -----------
# This script defines a set of mouse and human atlas labels that can be
# used to prune the neuroanatomical trees to the desired level of granularity

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.tree))
suppressPackageStartupMessages(library(optparse))

# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--outdir",
              type = "character",
              default = "data/",
              help = paste("Directory in which to save the .RData file",
                           "containing the label sets. [default %default")),
  make_option("--mousetree",
              type = "character",
              help = "Path to mouse tree .RData file."),
  make_option("--humantree",
              type = "character",
              help = "Path to human tree .RData file.")
)

args <- parse_args(OptionParser(option_list = option_list))

# Functions ------------------------------------------------------------------

working_dir <- getwd()

script_dir <- commandArgs() %>% 
  str_subset("--file=") %>% 
  str_remove("--file=") %>% 
  dirname()

path_tools <- file.path(working_dir, script_dir, "functions", "tree_tools.R")

souce(path_tools)

# Mouse labels 134 / Human labels 166 -------------------------------------------

#Import AMBA tree

load(args[["mousetree"]])

#Remove white matter and ventricles
pruneAnatTree(treeMouseExpr,
              nodes = c("fiber tracts", "ventricular systems"), 
              method = "AtNode")

#Remove bilateral regions
Prune(treeMouseExpr, 
      pruneFun = function(node){!str_detect(node$name, "left|right")})

labelsMouse_134 <- treeMouseExpr$Get("name", filterFun = isLeaf)
names(labelsMouse_134) <- NULL

#Import AHBA data tree
load(args[["humantree"]])

#Remove white matter and ventricles
pruneAnatTree(treeHumanExpr,
              nodes = c("white matter", "sulci & spaces"),
              method = "AtNode")

#Remove bilateral regions
Prune(treeHumanExpr,
      pruneFun = function(node){!str_detect(node$name, "left|right")})

labelsHuman_166 <- treeHumanExpr$Get("name", filterFun = isLeaf)
names(labelsHuman_166) <- NULL


# Mouse labels 67 / Human labels 88 -------------------------------------------
labelsMouse_67 <- c("Cortical subplate-other",
                    "Claustrum",
                    "Endopiriform nucleus",
                    "Olfactory areas-other",
                    "Postpiriform transition area",
                    "Piriform area",
                    "Taenia tecta",
                    "Cortical amygdalar area",
                    "Piriform-amygdalar area",
                    "Main olfactory bulb",
                    "Accessory olfactory bulb",
                    "Anterior olfactory nucleus",
                    "Subiculum",
                    "Entorhinal area",
                    "Field CA1",
                    "Field CA2",
                    "Field CA3",
                    "Dentate gyrus",
                    "Anterior cingulate area",
                    "Infralimbic area",
                    "Retrosplenial area",
                    "Prelimbic area",
                    "Primary auditory area",
                    "Dorsal auditory area",
                    "Ventral auditory area",
                    "Agranular insular area",
                    "Ectorhinal area",
                    "Primary motor area",
                    "Secondary motor area",
                    "Frontal pole, cerebral cortex",
                    "Orbital area",
                    "Posterior parietal association areas",
                    "Perirhinal area",
                    "Primary somatosensory area",
                    "Supplemental somatosensory area",
                    "Temporal association areas",
                    "Visual areas",
                    "Pallidum",
                    "Striatum ventral region",
                    "Lateral septal complex",
                    "Caudoputamen",
                    "Medial amygdalar nucleus",
                    "Inferior colliculus",
                    "Superior colliculus, sensory related",
                    "Midbrain raphe nuclei",
                    "Midbrain-other",
                    "Periaqueductal gray",
                    "Medulla",
                    "Pons",
                    "Hypothalamus",
                    "Thalamus",
                    "Lingula (I)",
                    "Central lobule",
                    "Culmen",
                    "Declive (VI)",
                    "Folium-tuber vermis (VII)",
                    "Pyramus (VIII)",
                    "Uvula (IX)",
                    "Nodulus (X)",
                    "Simple lobule",
                    "Crus 1",
                    "Crus 2",
                    "Paramedian lobule",
                    "Copula pyramidis",
                    "Flocculus",
                    "Paraflocculus",
                    "Cerebellar nuclei")

labelsHuman_88 <- c("anterior orbital gyrus",
                    "frontal operculum",
                    "frontal pole",
                    "gyrus rectus",
                    "inferior frontal gyrus",
                    "inferior rostral gyrus",
                    "lateral orbital gyrus",
                    "medial orbital gyrus",
                    "middle frontal gyrus",
                    "paracentral lobule, anterior part",
                    "paraterminal gyrus",
                    "parolfactory gyri",
                    "posterior orbital gyrus",
                    "precentral gyrus",
                    "superior frontal gyrus",
                    "superior rostral gyrus",
                    "insula",
                    "cingulate gyrus",
                    "dentate gyrus",
                    "CA1 field",
                    "CA2 field",
                    "CA3 field",
                    "CA4 field",
                    "subiculum",
                    "parahippocampal gyrus",
                    "piriform cortex",
                    "cuneus",
                    "inferior occipital gyrus",
                    "lingual gyrus",
                    "occipital pole",
                    "occipito-temporal gyrus",
                    "superior occipital gyrus",
                    "inferior parietal lobule",
                    "paracentral lobule, posterior part",
                    "postcentral gyrus",
                    "precuneus",
                    "superior parietal lobule",
                    "fusiform gyrus",
                    "Heschl's gyrus",
                    "inferior temporal gyrus",
                    "middle temporal gyrus",
                    "planum polare",
                    "planum temporale",
                    "superior temporal gyrus",
                    "temporal pole",
                    "transverse gyri",
                    "amygdala",
                    "septal nuclei",
                    "substantia innominata",
                    "globus pallidus",
                    "caudate nucleus",
                    "nucleus accumbens",
                    "putamen",
                    "claustrum",
                    "epithalamus",
                    "hypothalamus",
                    "subthalamus",
                    "thalamus",
                    "inferior colliculus",
                    "superior colliculus",
                    "midbrain tegmentum",
                    "pretectal region",
                    "vermal I-II",
                    "vermal III",
                    "vermal IV",
                    "vermal V",
                    "vermal VI",
                    "vermal VIIAf",
                    "vermal VIIAt",
                    "vermal VIIB",
                    "vermal VIIIA",
                    "vermal VIIIB",
                    "vermal IX",
                    "vermal X",
                    "III",
                    "IV",
                    "V",
                    "VI",
                    "crus I",
                    "crus II",
                    "VIIB",
                    "VIIIA",
                    "VIIIB",
                    "IX",
                    "X",
                    "cerebellar nuclei",
                    "pons",
                    "myelencephalon")

# Mouse labels 46 / Human labels 79 -------------------------------------------
labelsMouse_46 <- c("Cortical subplate-other",
                    "Claustrum",
                    "Endopiriform nucleus",
                    "Olfactory areas",
                    "Retrohippocampal region",
                    "Hippocampal region",
                    "Anterior cingulate area",
                    "Infralimbic area",
                    "Retrosplenial area",
                    "Prelimbic area",
                    "Auditory areas",
                    "Agranular insular area",
                    "Ectorhinal area",
                    "Somatomotor areas",
                    "Frontal pole, cerebral cortex",
                    "Orbital area",
                    "Posterior parietal association areas",
                    "Perirhinal area",
                    "Somatosensory areas",
                    "Temporal association areas",
                    "Visual areas",
                    "Pallidum",
                    "Striatum",
                    "Midbrain, sensory related",
                    "Midbrain, behavioral state related",
                    "Midbrain-other",
                    "Midbrain, motor related",
                    "Medulla",
                    "Pons",
                    "Hypothalamus",
                    "Thalamus",
                    "Lingula (I)",
                    "Central lobule",
                    "Culmen",
                    "Declive (VI)",
                    "Folium-tuber vermis (VII)",
                    "Pyramus (VIII)",
                    "Uvula (IX)",
                    "Nodulus (X)",
                    "Simple lobule",
                    "Ansiform lobule",
                    "Paramedian lobule",
                    "Copula pyramidis",
                    "Flocculus",
                    "Paraflocculus",
                    "Cerebellar nuclei")

labelsHuman_79 <- c("anterior orbital gyrus",
                    "frontal operculum",
                    "frontal pole",
                    "gyrus rectus",
                    "inferior frontal gyrus",
                    "inferior rostral gyrus",
                    "lateral orbital gyrus",
                    "medial orbital gyrus",
                    "middle frontal gyrus",
                    "paracentral lobule, anterior part",
                    "paraterminal gyrus",
                    "parolfactory gyri",
                    "posterior orbital gyrus",
                    "precentral gyrus",
                    "superior frontal gyrus",
                    "superior rostral gyrus",
                    "insula",
                    "cingulate gyrus",
                    "hippocampal formation",
                    "parahippocampal gyrus",
                    "piriform cortex",
                    "cuneus",
                    "inferior occipital gyrus",
                    "lingual gyrus",
                    "occipital pole",
                    "occipito-temporal gyrus",
                    "superior occipital gyrus",
                    "inferior parietal lobule",
                    "paracentral lobule, posterior part",
                    "postcentral gyrus",
                    "precuneus",
                    "superior parietal lobule",
                    "fusiform gyrus",
                    "Heschl's gyrus",
                    "inferior temporal gyrus",
                    "middle temporal gyrus",
                    "planum polare",
                    "planum temporale",
                    "superior temporal gyrus",
                    "temporal pole",
                    "transverse gyri",
                    "amygdala",
                    "basal forebrain",
                    "globus pallidus",
                    "striatum",
                    "claustrum",
                    "epithalamus",
                    "hypothalamus",
                    "subthalamus",
                    "thalamus",
                    "midbrain tectum",
                    "midbrain tegmentum",
                    "pretectal region",
                    "vermal I-II",
                    "vermal III",
                    "vermal IV",
                    "vermal V",
                    "vermal VI",
                    "vermal VIIAf",
                    "vermal VIIAt",
                    "vermal VIIB",
                    "vermal VIIIA",
                    "vermal VIIIB",
                    "vermal IX",
                    "vermal X",
                    "III",
                    "IV",
                    "V",
                    "VI",
                    "crus I",
                    "crus II",
                    "VIIB",
                    "VIIIA",
                    "VIIIB",
                    "IX",
                    "X",
                    "cerebellar nuclei",
                    "pons",
                    "myelencephalon")

# Mouse labels 28 / Human labels 56 -------------------------------------------
labelsMouse_28 <- c("Cortical subplate",
                    "Olfactory areas",
                    "Hippocampal formation",
                    "Anterior cingulate area",
                    "Infralimbic area",
                    "Retrosplenial area",
                    "Prelimbic area",
                    "Auditory areas",
                    "Agranular insular area",
                    "Ectorhinal area",
                    "Somatomotor areas",
                    "Frontal pole, cerebral cortex",
                    "Orbital area",
                    "Posterior parietal association areas",
                    "Perirhinal area",
                    "Somatosensory areas",
                    "Temporal association areas",
                    "Visual areas",
                    "Pallidum",
                    "Striatum",
                    "Midbrain",
                    "Medulla",
                    "Pons",
                    "Hypothalamus",
                    "Thalamus",
                    "Vermal regions",
                    "Hemispheric regions",
                    "Cerebellar nuclei")

labelsHuman_56 <- c("anterior orbital gyrus",
                    "frontal operculum",
                    "frontal pole",
                    "gyrus rectus",
                    "inferior frontal gyrus",
                    "inferior rostral gyrus",
                    "lateral orbital gyrus",
                    "medial orbital gyrus",
                    "middle frontal gyrus",
                    "paracentral lobule, anterior part",
                    "paraterminal gyrus",
                    "parolfactory gyri",
                    "posterior orbital gyrus",
                    "precentral gyrus",
                    "superior frontal gyrus",
                    "superior rostral gyrus",
                    "insula",
                    "cingulate gyrus",
                    "hippocampal formation",
                    "parahippocampal gyrus",
                    "piriform cortex",
                    "cuneus",
                    "inferior occipital gyrus",
                    "lingual gyrus",
                    "occipital pole",
                    "occipito-temporal gyrus",
                    "superior occipital gyrus",
                    "inferior parietal lobule",
                    "paracentral lobule, posterior part",
                    "postcentral gyrus",
                    "precuneus",
                    "superior parietal lobule",
                    "fusiform gyrus",
                    "Heschl's gyrus",
                    "inferior temporal gyrus",
                    "middle temporal gyrus",
                    "planum polare",
                    "planum temporale",
                    "superior temporal gyrus",
                    "temporal pole",
                    "transverse gyri",
                    "amygdala",
                    "basal forebrain",
                    "globus pallidus",
                    "striatum",
                    "claustrum",
                    "epithalamus",
                    "hypothalamus",
                    "subthalamus",
                    "thalamus",
                    "mesencephalon",
                    "vermis",
                    "cerebellar hemispheres",
                    "cerebellar nuclei",
                    "pons",
                    "myelencephalon")

# Mouse labels 11 / Human labels 16 -------------------------------------------
labelsMouse_11 <- c("Cortical subplate",
                    "Olfactory areas",
                    "Hippocampal formation",
                    "Isocortex",
                    "Cerebral nuclei",
                    "Midbrain",
                    "Medulla",
                    "Pons",
                    "Interbrain",
                    "Cerebellar cortex",
                    "Cerebellar nuclei")

labelsHuman_16 <- c("frontal lobe",
                    "insula",
                    "limbic lobe",
                    "occipital lobe",
                    "parietal lobe",
                    "temporal lobe",
                    "amygdala",
                    "basal forebrain",
                    "basal ganglia",
                    "claustrum",
                    "diencephalon",
                    "mesencephalon",
                    "cerebellar cortex",
                    "cerebellar nuclei",
                    "pons",
                    "myelencephalon")


# Mouse labels 5 / Human labels 5 ---------------------------------------------
labelsMouse_5 <- c("Cerebrum", 
                   "Interbrain",
                   "Midbrain",
                   "Hindbrain",
                   "Cerebellum")

labelsHuman_5 <- c("telencephalon",
                   "diencephalon",
                   "mesencephalon",
                   "metencephalon",
                   "myelencephalon")


# Mouse labels 11 / Human labels 16 re-ordered -------------------------------------------
labelsMouse_11_reordered <- c("Cortical subplate",
                              "Olfactory areas",
                              "Hippocampal formation",
                              "Isocortex",
                              "Cerebral nuclei",
                              "Interbrain",
                              "Midbrain",
                              "Pons",
                              "Medulla",
                              "Cerebellar cortex",
                              "Cerebellar nuclei")

labelsHuman_16_reordered <- c("claustrum",
                              "limbic lobe",
                              "frontal lobe",
                              "insula",
                              "occipital lobe",
                              "parietal lobe",
                              "temporal lobe",
                              "amygdala",
                              "basal ganglia",
                              "basal forebrain",
                              "diencephalon",
                              "mesencephalon",
                              "pons",
                              "myelencephalon",
                              "cerebellar cortex",
                              "cerebellar nuclei")

#Create a tree with 67 leaf nodes
treeMouse_67 <- Clone(treeMouseExpr)
pruneAnatTree(treeMouse_67,
              nodes = labelsMouse_67,
              method = "BelowNode")

#Order the 67 labels according to the ordered 11 labels
labelsMouse_67_reordered <- character(0)
for(lab in labelsMouse_11_reordered){
  labelsMouse_67_reordered <- c(labelsMouse_67_reordered, FindNode(treeMouse_67, lab)$Get("name", filterFun = isLeaf))
}
names(labelsMouse_67_reordered) <- NULL

#Create a tree with 88 leaf nodes
treeHuman_88 <- Clone(treeHumanExpr)
pruneAnatTree(treeHuman_88,
              nodes = labelsHuman_88,
              method = "BelowNode")

#Order the 88 labels according to the ordered 16 labels
labelsHuman_88_reordered <- character(0)
for(lab in labelsHuman_16_reordered){
  labelsHuman_88_reordered <- c(labelsHuman_88_reordered, FindNode(treeHuman_88, lab)$Get("name", filterFun = isLeaf))
}
names(labelsHuman_88_reordered) <- NULL

#Move around limbic node regions to better match the mouse labels
ind <- c(1, 3:8, 2, 9:length(labelsHuman_88_reordered))
labelsHuman_88_reordered <- labelsHuman_88_reordered[ind]


# Label list definitions ------------------------------------------------------

listLabelsMouse <- list(Region5 = labelsMouse_5,
                        Region11 = labelsMouse_11,
                        Region28 = labelsMouse_28,
                        Region46 = labelsMouse_46,
                        Region67 = labelsMouse_67,
                        Region134 = labelsMouse_134)

listLabelsMouseReordered <- list(Region11_reordered = labelsMouse_11_reordered,
                                 Region67_reordered = labelsMouse_67_reordered)

listLabelsHuman = list(Region5 = labelsHuman_5,
                       Region16 = labelsHuman_16,
                       Region56 = labelsHuman_56,
                       Region79 = labelsHuman_79,
                       Region88 = labelsHuman_88,
                       Region166 = labelsHuman_166)

listLabelsHumanReordered <- list(Region16_reordered = labelsHuman_16_reordered,
                                 Region88_reordered = labelsHuman_88_reordered)

fileout <- file.path(args[["outdir"]], "TreeLabels.RData")
fileout_reordered <- file.path(args[["outdir"]], "TreeLabelsReordered.RData")

save(listLabelsMouse,
     listLabelsHuman,
     file = fileout)

save(listLabelsMouseReordered,
     listLabelsHumanReordered,
     file = fileout_reordered)