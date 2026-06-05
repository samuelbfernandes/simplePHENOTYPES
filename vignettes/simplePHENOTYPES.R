## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#",
  root.dir = getwd()
)

## ----load data, eval = FALSE--------------------------------------------------
#  library(simplePHENOTYPES)
#  data("SNP55K_maize282_maf04")
#  SNP55K_maize282_maf04[1:8, 1:10]

## ----ST, eval = FALSE---------------------------------------------------------
#  create_phenotypes(
#    geno_obj = SNP55K_maize282_maf04,
#    add_QTN_num = 3,
#    add_effect = 0.2,
#    big_add_QTN_effect = 0.9,
#    rep = 10,
#    h2 = 0.7,
#    model = "A",
#    home_dir = tempdir())

## ----MT P, results = "hide", eval = FALSE-------------------------------------
#   test1 <-  create_phenotypes(
#      geno_obj = SNP55K_maize282_maf04,
#      add_QTN_num = 3,
#      dom_QTN_num = 4,
#      big_add_QTN_effect = c(0.3, 0.3, 0.3),
#      h2 = c(0.2, 0.4, 0.4),
#      add_effect = c(0.04,0.2,0.1),
#      dom_effect = c(0.04,0.2,0.1),
#      ntraits = 3,
#      rep = 10,
#      vary_QTN = FALSE,
#      output_format = "multi-file",
#      architecture = "pleiotropic",
#      output_dir = "Results_Pleiotropic",
#      to_r = TRUE,
#      seed = 10,
#      model = "AD",
#      sim_method = "geometric",
#    home_dir = tempdir()
#    )

## ----MT P2, results = "hide", eval = FALSE------------------------------------
#   custom_geometric_a <- list(trait_1 = c(0.04, 0.0016),
#                           trait_2 = c(0.2, 0.04),
#                           trait_3 = c(0.1, 0.01))
#   custom_geometric_d <- list(trait_1 = c(0.04, 0.0016, 6.4e-05, 2.56e-06),
#                           trait_2 = c(0.2, 0.04, 0.008, 0.0016),
#                           trait_3 = c(0.1, 0.01, 0.001, 1e-04))
#  
#   test2 <-  create_phenotypes(
#     geno_obj = SNP55K_maize282_maf04,
#     add_QTN_num = 3,
#     dom_QTN_num = 4,
#     big_add_QTN_effect = c(0.3, 0.3, 0.3),
#     h2 = c(0.2,0.4, 0.4),
#     add_effect = custom_geometric_a,
#     dom_effect = custom_geometric_d,
#     ntraits = 3,
#     rep = 10,
#     vary_QTN = FALSE,
#     output_format = "multi-file",
#     architecture = "pleiotropic",
#     output_dir = "Results_Pleiotropic",
#     to_r = T,
#     sim_method = "custom",
#     seed = 10,
#     model = "AD",
#    home_dir = tempdir()
#   )
#  
#   all.equal(test1, test2)

## ----MT PP, results = "hide", eval = FALSE------------------------------------
#  cor_matrix <- matrix(c(   1, 0.3, -0.9,
#                          0.3,   1,  -0.5,
#                         -0.9, -0.5,    1 ), 3)
#  
#  sim_results <- create_phenotypes(
#    geno_obj = SNP55K_maize282_maf04,
#    ntraits = 3,
#    pleio_a = 3,
#    pleio_e = 2,
#    same_add_dom_QTN = TRUE,
#    degree_of_dom = 0.5,
#    trait_spec_a_QTN_num = c(4, 10, 1),
#    trait_spec_e_QTN_num = c(3, 2, 5),
#    h2 = c(0.2, 0.4, 0.8),
#    add_effect = c(0.5, 0.33, 0.2),
#    epi_effect = c(0.3, 0.3, 0.3),
#    epi_interaction = 2,
#    cor = cor_matrix,
#    rep = 20,
#    output_dir = "Results_Partially",
#    output_format = "long",
#    architecture = "partially",
#    out_geno = "numeric",
#    to_r = TRUE,
#    model = "AE",
#    home_dir = tempdir()
#  )

## ----MT LD, results = "hide", eval = FALSE------------------------------------
#  create_phenotypes(
#    geno_obj = SNP55K_maize282_maf04,
#    add_QTN_num = 3,
#    h2 = c(0.2, 0.4),
#    add_effect = c(0.02, 0.05),
#    rep = 5,
#    seed = 200,
#    output_format = "wide",
#    architecture = "LD",
#    output_dir = "Results_LD",
#    out_geno = "BED",
#    remove_QTN = TRUE,
#    ld_max =0.8,
#    ld_min =0.2,
#    model = "A",
#    ld_method = "composite",
#    type_of_ld = "indirect",
#    home_dir = tempdir()
#  )

## ----MT PP E, results = "hide", eval = FALSE----------------------------------
#  residual <- matrix(c(1, 0.1,-0.2,
#                       0.1, 1,-0.1,-0.2,-0.1, 1), 3)
#  heritability <- matrix(c(0.2, 0.4, 0.8,
#                           0.6, 0.7, 0.2), 2)
#  create_phenotypes(
#    geno_obj = SNP55K_maize282_maf04,
#    pleio_a = 3,
#    pleio_e = 2,
#    same_add_dom_QTN = TRUE,
#    degree_of_dom = 1,
#    trait_spec_a_QTN_num = c(4, 10, 1),
#    trait_spec_e_QTN_num = c(2, 1, 5),
#    epi_effect = c(0.01, 0.4, 0.2),
#    add_effect = c(0.3, 0.2, 0.5),
#    h2 = heritability,
#    ntraits = 3,
#    rep = 5,
#    vary_QTN = TRUE,
#    warning_file_saver = FALSE,
#    output_dir = "Results_Partially_ADE",
#    output_format = "gemma",
#    architecture = "partially",
#    model = "ADE",
#    QTN_variance = TRUE,
#    remove_QTN = TRUE,
#    home_dir = tempdir(),
#    constraints = list(
#      maf_above = 0.3,
#      maf_below = 0.44,
#      hets = "include"
#    ),
#    cor_res = residual
#  )

## ----qtn_list, results = "hide", eval = FALSE---------------------------------
#  QTN_list <- list()
#  QTN_list$add[[1]] <- c("ss196523212")
#  QTN_list$dom[[1]] <- c("ss196510214", "ss196472187")
#  QTN_list$epi[[1]] <- c("ss196530605", "ss196475446")
#  create_phenotypes(
#    geno_obj = SNP55K_maize282_maf04,
#    add_QTN_num = 1,
#    dom_QTN_num = 2,
#    epi_QTN_num = 1,
#    epi_interaction = 2,
#    h2 = c(0.92, 0.4) ,
#    add_effect = c(0.90, 0.2),
#    dom_effect = c(0.01, 0.3),
#    epi_effect = c(-0.3, 0.7),
#    ntraits = 2,
#    QTN_list = QTN_list,
#    rep = 1,
#    output_format = "gemma",
#    out_geno = "BED",
#    output_dir = "output_data",
#    model = "ADE",
#    home_dir = getwd()
#  )
#  

## ----example, results = "hide", eval = FALSE----------------------------------
#  create_phenotypes(
#    geno_path = "PATH/TO/FILE",
#    prefix = "WGS_chrm_",
#    add_QTN_num = 3,
#    h2 = 0.2,
#    add_effect = 0.02,
#    rep = 5,
#    seed = 200,
#    output_format = "gemma",
#    output_dir = "Results",
#    model = "ADE",
#    home_dir = tempdir()
#  )

