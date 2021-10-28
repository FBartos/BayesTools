context("Prior informed function")


test_that("Informed prior distributions match the specification", {

  ### psychology
  # Oosterwijk
  p1 <- prior_informed("Oosterwijk")
  expect_equal(p1$distribution, "t")
  expect_equal(p1$parameters, list("location" = .35, "scale" = .102, "df" = 3))
  expect_equal(p1$truncation, list("lower" = 0, "upper" = Inf))

  # van Erp
  p2 <- prior_informed("van Erp")
  expect_equal(p2$distribution, "invgamma")
  expect_equal(p2$parameters, list("shape" = 1, "scale" = .15))
  expect_equal(p2$truncation, list("lower" = 0, "upper" = Inf))

  ### medicine
  # table from the manuscript - with slight formatting edits for computer readability
  medicine_table <-
  "Acute Respiratory Infections , 6 , 104 , mathcal{T}(,0, 0.38, 5,) , text{Inv-Gamma}(,1.73, 0.46,)
  ,Airways , 46 , 815 , mathcal{T}(,0, 0.38, 6,) , text{Inv-Gamma}(,2.02, 0.28,)
  ,Anaesthesia , 44 , 661 , mathcal{T}(,0, 0.55, 4,) , text{Inv-Gamma}(,1.62, 0.64,)
  ,Back and Neck , 13 , 278 , mathcal{T}(,0, 0.37, 5,) , text{Inv-Gamma}(,1.75, 0.57,)
  ,Bone Joint and Muscle Trauma , 32 , 1221 , mathcal{T}(,0, 0.40, 5,) , text{Inv-Gamma}(,1.52, 0.28,)
  ,Colorectal , 13 , 372 , mathcal{T}(,0, 0.51, 5,) , text{Inv-Gamma}(,1.64, 0.56,)
  ,Common Mental Disorders , 17 , 264 , mathcal{T}(,0, 0.55, 5,) , text{Inv-Gamma}(,1.62, 0.45,)
  ,Consumers and Communication , 6 , 72 , mathcal{T}(,0, 0.40, 5,) , text{Inv-Gamma}(,1.56, 0.14,)
  ,Cystic Fibrosis and Genetic Disorders , 1 , 12 , mathcal{T}(,0, 0.47, 5,) , text{Inv-Gamma}(,1.70, 0.45,)
  ,Dementia and Cognitive Improvement , 9 , 197 , mathcal{T}(,0, 0.45, 5,) , text{Inv-Gamma}(,1.71, 0.44,)
  ,Developmental Psychosocial and Learning Problems , 20 , 407 , mathcal{T}(,0, 0.18, 5,) , text{Inv-Gamma}(,1.43, 0.12,)
  ,Drugs and Alcohol , 8 , 170 , mathcal{T}(,0, 0.33, 5,) , text{Inv-Gamma}(,1.89, 0.28,)
  ,Effective Practice and Organisation of Care , 10 , 204 , mathcal{T}(,0, 0.39, 5,) , text{Inv-Gamma}(,1.71, 0.35,)
  ,Emergency and Critical Care , 9 , 214 , mathcal{T}(,0, 0.39, 5,) , text{Inv-Gamma}(,1.62, 0.29,)
  ,ENT , 17 , 273 , mathcal{T}(,0, 0.43, 5,) , text{Inv-Gamma}(,1.85, 0.48,)
  ,Eyes and Vision , 14 , 347 , mathcal{T}(,0, 0.40, 6,) , text{Inv-Gamma}(,1.86, 0.41,)
  ,Gynaecological Neuro-oncology and Orphan Cancer , 1 , 10 , mathcal{T}(,0, 0.45, 5,) , text{Inv-Gamma}(,1.67, 0.46,)
  ,Gynaecology and Fertility , 14 , 253 , mathcal{T}(,0, 0.38, 5,) , text{Inv-Gamma}(,1.78, 0.46,)
  ,Heart , 88 , 2112 , mathcal{T}(,0, 0.42, 5,) , text{Inv-Gamma}(,1.83, 0.47,)
  ,Hepato-Biliary , 34 , 1103 , mathcal{T}(,0, 0.60, 4,) , text{Inv-Gamma}(,1.56, 0.58,)
  ,HIV/AIDS , 2 , 23 , mathcal{T}(,0, 0.43, 5,) , text{Inv-Gamma}(,1.73, 0.44,)
  ,Hypertension , 27 , 524 , mathcal{T}(,0, 0.48, 3,) , text{Inv-Gamma}(,2.01, 0.38,)
  ,Incontinence , 17 , 219 , mathcal{T}(,0, 0.33, 6,) , text{Inv-Gamma}(,1.64, 0.36,)
  ,Infectious Diseases , 8 , 150 , mathcal{T}(,0, 0.59, 2,) , text{Inv-Gamma}(,1.28, 0.44,)
  ,Inflammatory Bowel Disease , 1 , 12 , mathcal{T}(,0, 0.40, 5,) , text{Inv-Gamma}(,1.76, 0.39,)
  ,Injuries , 3 , 54 , mathcal{T}(,0, 0.35, 5,) , text{Inv-Gamma}(,1.80, 0.34,)
  ,Kidney and Transplant , 39 , 767 , mathcal{T}(,0, 0.54, 5,) , text{Inv-Gamma}(,1.72, 0.53,)
  ,Metabolic and Endocrine Disorders , 25 , 503 , mathcal{T}(,0, 0.43, 5,) , text{Inv-Gamma}(,1.71, 0.37,)
  ,Methodology , 5 , 106 , mathcal{T}(,0, 0.49, 5,) , text{Inv-Gamma}(,1.72, 0.51,)
  ,Movement Disorders , 5 , 70 , mathcal{T}(,0, 0.42, 5,) , text{Inv-Gamma}(,1.88, 0.33,)
  ,Musculoskeletal , 32 , 778 , mathcal{T}(,0, 0.45, 6,) , text{Inv-Gamma}(,1.87, 0.38,)
  ,Neonatal , 11 , 259 , mathcal{T}(,0, 0.42, 5,) , text{Inv-Gamma}(,1.68, 0.38,)
  ,Oral Health , 10 , 236 , mathcal{T}(,0, 0.51, 5,) , text{Inv-Gamma}(,1.79, 0.28,)
  ,Pain Palliative and Supportive Care , 16 , 283 , mathcal{T}(,0, 0.43, 5,) , text{Inv-Gamma}(,1.69, 0.42,)
  ,Pregnancy and Childbirth , 32 , 539 , mathcal{T}(,0, 0.33, 5,) , text{Inv-Gamma}(,1.86, 0.32,)
  ,Public Health , 2 , 22 , mathcal{T}(,0, 0.33, 5,) , text{Inv-Gamma}(,1.76, 0.23,)
  ,Schizophrenia , 21 , 436 , mathcal{T}(,0, 0.29, 4,) , text{Inv-Gamma}(,1.60, 0.27,)
  ,Sexually Transmitted Infections , 9 , 113 , mathcal{T}(,0, 0.42, 5,) , text{Inv-Gamma}(,1.70, 0.59,)
  ,Skin , 6 , 85 , mathcal{T}(,0, 0.48, 5,) , text{Inv-Gamma}(,1.64, 0.51,)
  ,Stroke , 21 , 357 , mathcal{T}(,0, 0.48, 5,) , text{Inv-Gamma}(,1.71, 0.40,)
  ,Tobacco Addiction , 4 , 44 , mathcal{T}(,0, 0.44, 4,) , text{Inv-Gamma}(,1.73, 0.42,)
  ,Upper GI and Pancreatic Diseases , 1 , 12 , mathcal{T}(,0, 0.45, 5,) , text{Inv-Gamma}(,1.76, 0.38,)
  ,Urology , 2 , 33 , mathcal{T}(,0, 0.44, 5,) , text{Inv-Gamma}(,1.73, 0.45,)
  ,Vascular , 3 , 35 , mathcal{T}(,0, 0.46, 5,) , text{Inv-Gamma}(,1.66, 0.50,)
  ,Work , 2 , 24 , mathcal{T}(,0, 0.42, 5,) , text{Inv-Gamma}(,1.76, 0.39,)
  ,Wounds , 7 , 103 , mathcal{T}(,0, 0.56, 5,) , text{Inv-Gamma}(,1.54, 0.41,)
  ,Cochrane , 713 , 14876 , mathcal{T}(,0, 0.43, 5,) , text{Inv-Gamma}(,1.71, 0.40,)"
  medicine_table <- strsplit(medicine_table, ",")[[1]]
  medicine_table <- matrix(medicine_table, ncol = 12, byrow = T)

  for(i in 1:nrow(medicine_table)){
    # test prior for the effect
    p1 <- prior_informed(medicine_table[i,1], parameter = "effect", type = "smd")
    expect_equal(p1$distribution, "t")
    expect_equal(p1$parameters, list("location" = as.numeric(medicine_table[i,5]), "scale" = as.numeric(medicine_table[i,6]), "df" = as.numeric(medicine_table[i,7])))
    expect_equal(p1$truncation, list("lower" = -Inf, "upper" = Inf))

    # test prior for the heterogeneity
    p2 <- prior_informed(medicine_table[i,1], parameter = "heterogeneity", type = "smd")
    expect_equal(p2$distribution, "invgamma")
    expect_equal(p2$parameters, list("shape" = as.numeric(medicine_table[i,10]), "scale" = as.numeric(medicine_table[i,11])))
    expect_equal(p2$truncation, list("lower" = 0, "upper" = Inf))
  }

  ### other
  expect_error(prior_informed("random"))
  expect_error(prior_informed("Cochrane", parameter = "random"))
  expect_error(prior_informed("Cochrane", parameter = "effect", type = "logOR"))

})
