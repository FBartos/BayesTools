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



  ### dichotomous and time to event
  paper_logOR <-  "
  Acute Respiratory Infections                              & Student-t(0, 0.48, 3) & Inv-Gamma(1.67, 0.45)  &
  Airways                                                   & Student-t(0, 0.37, 2) & Inv-Gamma(1.35, 0.27)  &
  Anaesthesia                                               & Student-t(0, 0.79, 6) & Inv-Gamma(2.12, 0.86)  &
  Back and Neck                                             & Student-t(0, 0.62, 4) & Inv-Gamma(1.84, 0.68)  &
  Bone, Joint and Muscle Trauma                             & Student-t(0, 0.44, 2) & Inv-Gamma(1.01, 0.36)  &
  Breast Cancer                                             & Student-t(0, 0.39, 3) & Inv-Gamma(1.59, 0.48)  &
  Childhood Cancer                                          & Student-t(0, 0.47, 4) & ---                    &
  Colorectal                                                & Student-t(0, 0.65, 4) & Inv-Gamma(1.71, 0.60)  &
  Common Mental Disorders                                   & Student-t(0, 0.54, 4) & Inv-Gamma(1.63, 0.45)  &
  Consumers and Communication                               & Student-t(0, 0.48, 4) & Inv-Gamma(2.37, 0.86)  &
  Cystic Fibrosis and Genetic Disorders                     & Student-t(0, 0.40, 3) & Inv-Gamma(1.82, 0.58)  &
  Dementia and Cognitive Improvement                        & Student-t(0, 0.49, 4) & Inv-Gamma(1.28, 0.25)  &
  Developmental, Psychosocial and Learning Problems         & Student-t(0, 0.83, 5) & Inv-Gamma(1.83, 0.82)  &
  Drugs and Alcohol                                         & Student-t(0, 0.44, 4) & Inv-Gamma(1.59, 0.42)  &
  Effective Practice and Organisation of Care               & Student-t(0, 0.51, 4) & Inv-Gamma(1.90, 0.68)  &
  Emergency and Critical Care                               & Student-t(0, 0.35, 3) & Inv-Gamma(1.46, 0.34)  &
  ENT                                                       & Student-t(0, 0.81, 4) & Inv-Gamma(1.73, 0.71)  &
  Epilepsy                                                  & Student-t(0, 0.88, 6) & Inv-Gamma(1.71, 0.43)  &
  Eyes and Vision                                           & Student-t(0, 0.77, 5) & Inv-Gamma(2.09, 0.94)  &
  Fertility Regulation                                      & Student-t(0, 0.46, 5) & Inv-Gamma(2.21, 0.71)  &
  Gut                                                       & Student-t(0, 0.63, 5) & Inv-Gamma(1.94, 0.62)  &
  Gynaecological, Neuro-oncology and Orphan Cancer          & Student-t(0, 0.53, 4) & Inv-Gamma(1.80, 0.56)  &
  Gynaecology and Fertility                                 & Student-t(0, 0.40, 2) & Inv-Gamma(1.24, 0.28)  &
  Haematology                                               & Student-t(0, 0.57, 4) & Inv-Gamma(2.91, 0.66)  &
  Heart                                                     & Student-t(0, 0.20, 2) & Inv-Gamma(1.64, 0.29)  &
  Heart; Vascular                                           & Student-t(0, 0.95, 4) & Inv-Gamma(1.64, 0.83)  &
  Hepato-Biliary                                            & Student-t(0, 0.43, 3) & Inv-Gamma(1.58, 0.40)  &
  HIV/AIDS                                                  & Student-t(0, 0.32, 4) & Inv-Gamma(1.76, 0.36)  &
  Hypertension                                              & Student-t(0, 0.28, 5) & Inv-Gamma(1.25, 0.10)  &
  Incontinence                                              & Student-t(0, 0.75, 3) & Inv-Gamma(2.07, 1.09)  &
  Infectious Diseases                                       & Student-t(0, 0.66, 3) & Inv-Gamma(2.08, 0.86)  &
  Injuries                                                  & Student-t(0, 0.60, 4) & Inv-Gamma(1.52, 0.49)  &
  Kidney and Transplant                                     & Student-t(0, 0.53, 4) & Inv-Gamma(1.68, 0.44)  &
  Lung Cancer                                               & Student-t(0, 0.61, 5) & Inv-Gamma(2.04, 0.68)  &
  Metabolic and Endocrine Disorders                         & Student-t(0, 0.29, 2) & Inv-Gamma(0.92, 0.11)  &
  Methodology                                               & Student-t(0, 0.60, 5) & Inv-Gamma(2.04, 0.50)  &
  Movement Disorders                                        & Student-t(0, 0.73, 5) & Inv-Gamma(2.14, 0.64)  &
  Multiple Sclerosis and Rare Diseases of the CNS           & Student-t(0, 0.76, 4) & Inv-Gamma(2.09, 0.71)  &
  Musculoskeletal                                           & Student-t(0, 0.59, 4) & Inv-Gamma(1.76, 0.59)  &
  Neonatal                                                  & Student-t(0, 0.29, 3) & Inv-Gamma(1.80, 0.42)  &
  Neuromuscular                                             & Student-t(0, 0.70, 5) & Inv-Gamma(1.74, 0.49)  &
  Oral Health                                               & Student-t(0, 1.13, 4) & Inv-Gamma(1.85, 0.70)  &
  Pain, Palliative and Supportive Care                      & Student-t(0, 1.00, 7) & Inv-Gamma(1.50, 0.40)  &
  Pregnancy and Childbirth                                  & Student-t(0, 0.38, 2) & Inv-Gamma(1.73, 0.47)  &
  Schizophrenia                                             & Student-t(0, 0.58, 4) & Inv-Gamma(1.92, 0.69)  &
  Sexually Transmitted Infections                           & Student-t(0, 0.61, 2) & Inv-Gamma(1.72, 0.47)  &
  Skin                                                      & Student-t(0, 0.81, 2) & Inv-Gamma(1.64, 0.49)  &
  Stroke                                                    & Student-t(0, 0.22, 2) & Inv-Gamma(1.36, 0.21)  &
  Tobacco Addiction                                         & Student-t(0, 0.49, 5) & Inv-Gamma(2.51, 0.63)  &
  Urology                                                   & Student-t(0, 0.82, 5) & Inv-Gamma(1.72, 0.50)  &
  Vascular                                                  & Student-t(0, 0.68, 5) & Inv-Gamma(1.70, 0.45)  &
  Work                                                      & Student-t(0, 0.59, 4) & Inv-Gamma(1.79, 0.73)  &
  Wounds                                                    & Student-t(0, 0.60, 4) & Inv-Gamma(2.16, 0.86)  &
  Cochrane                                                  & Student-t(0, 0.58, 4) & Inv-Gamma(1.77, 0.55) "

  paper_logRR <-  "
  Acute Respiratory Infections                              & Student-t(0, 0.27, 3) & Inv-Gamma(1.58, 0.25) &
  Airways                                                   & Student-t(0, 0.25, 3) & Inv-Gamma(1.36, 0.15) &
  Anaesthesia                                               & Student-t(0, 0.45, 5) & Inv-Gamma(1.21, 0.21) &
  Back and Neck                                             & Student-t(0, 0.34, 3) & Inv-Gamma(1.52, 0.22) &
  Bone, Joint and Muscle Trauma                             & Student-t(0, 0.26, 2) & Inv-Gamma(0.97, 0.14) &
  Breast Cancer                                             & Student-t(0, 0.21, 3) & Inv-Gamma(1.49, 0.24) &
  Childhood Cancer                                          & Student-t(0, 0.28, 3) & ---                   &
  Colorectal                                                & Student-t(0, 0.38, 4) & Inv-Gamma(1.42, 0.27) &
  Common Mental Disorders                                   & Student-t(0, 0.29, 3) & Inv-Gamma(1.43, 0.20) &
  Consumers and Communication                               & Student-t(0, 0.13, 2) & Inv-Gamma(1.48, 0.14) &
  Cystic Fibrosis and Genetic Disorders                     & Student-t(0, 0.21, 2) & Inv-Gamma(1.71, 0.24) &
  Dementia and Cognitive Improvement                        & Student-t(0, 0.29, 3) & Inv-Gamma(1.53, 0.19) &
  Developmental, Psychosocial and Learning Problems         & Student-t(0, 0.39, 3) & Inv-Gamma(0.97, 0.11) &
  Drugs and Alcohol                                         & Student-t(0, 0.22, 4) & Inv-Gamma(1.58, 0.20) &
  Effective Practice and Organisation of Care               & Student-t(0, 0.30, 4) & Inv-Gamma(1.40, 0.24) &
  Emergency and Critical Care                               & Student-t(0, 0.19, 2) & Inv-Gamma(1.65, 0.25) &
  ENT                                                       & Student-t(0, 0.40, 4) & Inv-Gamma(1.50, 0.24) &
  Epilepsy                                                  & Student-t(0, 0.58, 5) & Inv-Gamma(1.74, 0.27) &
  Eyes and Vision                                           & Student-t(0, 0.34, 4) & Inv-Gamma(1.60, 0.26) &
  Fertility Regulation                                      & Student-t(0, 0.34, 4) & Inv-Gamma(1.59, 0.31) &
  Gut                                                       & Student-t(0, 0.40, 4) & Inv-Gamma(1.69, 0.29) &
  Gynaecological, Neuro-oncology and Orphan Cancer          & Student-t(0, 0.35, 4) & Inv-Gamma(1.43, 0.25) &
  Gynaecology and Fertility                                 & Student-t(0, 0.25, 2) & Inv-Gamma(1.29, 0.16) &
  Haematology                                               & Student-t(0, 0.31, 3) & Inv-Gamma(1.22, 0.11) &
  Heart                                                     & Student-t(0, 0.17, 3) & Inv-Gamma(1.96, 0.24) &
  Heart; Vascular                                           & Student-t(0, 0.39, 3) & Inv-Gamma(1.62, 0.25) &
  Hepato-Biliary                                            & Student-t(0, 0.23, 3) & Inv-Gamma(1.44, 0.19) &
  HIV/AIDS                                                  & Student-t(0, 0.21, 3) & Inv-Gamma(1.66, 0.23) &
  Hypertension                                              & Student-t(0, 0.25, 4) & Inv-Gamma(1.48, 0.11) &
  Incontinence                                              & Student-t(0, 0.39, 3) & Inv-Gamma(1.53, 0.33) &
  Infectious Diseases                                       & Student-t(0, 0.41, 3) & Inv-Gamma(1.30, 0.27) &
  Injuries                                                  & Student-t(0, 0.27, 4) & Inv-Gamma(1.46, 0.19) &
  Kidney and Transplant                                     & Student-t(0, 0.35, 4) & Inv-Gamma(1.51, 0.23) &
  Lung Cancer                                               & Student-t(0, 0.38, 4) & Inv-Gamma(1.49, 0.26) &
  Metabolic and Endocrine Disorders                         & Student-t(0, 0.13, 1) & Inv-Gamma(1.55, 0.19) &
  Methodology                                               & Student-t(0, 0.27, 5) & Inv-Gamma(1.74, 0.23) &
  Movement Disorders                                        & Student-t(0, 0.46, 4) & Inv-Gamma(1.54, 0.24) &
  Multiple Sclerosis and Rare Diseases of the CNS           & Student-t(0, 0.44, 3) & Inv-Gamma(1.58, 0.35) &
  Musculoskeletal                                           & Student-t(0, 0.33, 3) & Inv-Gamma(1.54, 0.31) &
  Neonatal                                                  & Student-t(0, 0.18, 3) & Inv-Gamma(1.89, 0.30) &
  Neuromuscular                                             & Student-t(0, 0.41, 4) & Inv-Gamma(1.44, 0.15) &
  Oral Health                                               & Student-t(0, 0.55, 4) & Inv-Gamma(1.50, 0.37) &
  Pain, Palliative and Supportive Care                      & Student-t(0, 0.59, 5) & Inv-Gamma(1.98, 0.39) &
  Pregnancy and Childbirth                                  & Student-t(0, 0.24, 2) & Inv-Gamma(1.40, 0.21) &
  Schizophrenia                                             & Student-t(0, 0.33, 3) & Inv-Gamma(1.47, 0.26) &
  Sexually Transmitted Infections                           & Student-t(0, 0.39, 3) & Inv-Gamma(1.43, 0.21) &
  Skin                                                      & Student-t(0, 0.41, 1) & Inv-Gamma(1.49, 0.23) &
  Stroke                                                    & Student-t(0, 0.13, 2) & Inv-Gamma(1.52, 0.15) &
  Tobacco Addiction                                         & Student-t(0, 0.36, 5) & Inv-Gamma(2.02, 0.36) &
  Urology                                                   & Student-t(0, 0.44, 3) & Inv-Gamma(1.37, 0.18) &
  Vascular                                                  & Student-t(0, 0.43, 5) & Inv-Gamma(1.26, 0.19) &
  Work                                                      & Student-t(0, 0.23, 3) & Inv-Gamma(1.50, 0.26) &
  Wounds                                                    & Student-t(0, 0.37, 4) & Inv-Gamma(1.63, 0.38) &
  Cochrane                                                  & Student-t(0, 0.32, 3) & Inv-Gamma(1.51, 0.23) "

  paper_RD <-  "
  Acute Respiratory Infections                              & Student-t(0, 0.01, 1) & Normal(0, 0.10) &
  Airways                                                   & Student-t(0, 0.01, 1) & Normal(0, 0.10) &
  Anaesthesia                                               & Student-t(0, 0.02, 1) & Normal(0, 0.11) &
  Back and Neck                                             & Student-t(0, 0.06, 2) & Normal(0, 0.10) &
  Bone, Joint and Muscle Trauma                             & Student-t(0, 0.01, 1) & Normal(0, 0.12) &
  Breast Cancer                                             & Student-t(0, 0.05, 2) & Normal(0, 0.12) &
  Childhood Cancer                                          & Student-t(0, 0.02, 1) & ---             &
  Colorectal                                                & Student-t(0, 0.01, 1) & Normal(0, 0.10) &
  Common Mental Disorders                                   & Student-t(0, 0.05, 1) & Normal(0, 0.11) &
  Consumers and Communication                               & Student-t(0, 0.06, 2) & Normal(0, 0.10) &
  Cystic Fibrosis and Genetic Disorders                     & Student-t(0, 0.01, 1) & Normal(0, 0.09) &
  Dementia and Cognitive Improvement                        & Student-t(0, 0.02, 1) & Normal(0, 0.09) &
  Developmental, Psychosocial and Learning Problems         & Student-t(0, 0.08, 2) & Normal(0, 0.13) &
  Drugs and Alcohol                                         & Student-t(0, 0.02, 1) & Normal(0, 0.10) &
  Effective Practice and Organisation of Care               & Student-t(0, 0.06, 2) & Normal(0, 0.11) &
  Emergency and Critical Care                               & Student-t(0, 0.02, 1) & Normal(0, 0.07) &
  ENT                                                       & Student-t(0, 0.04, 1) & Normal(0, 0.12) &
  Epilepsy                                                  & Student-t(0, 0.05, 2) & Normal(0, 0.09) &
  Eyes and Vision                                           & Student-t(0, 0.04, 1) & Normal(0, 0.11) &
  Fertility Regulation                                      & Student-t(0, 0.01, 1) & Normal(0, 0.08) &
  Gut                                                       & Student-t(0, 0.05, 2) & Normal(0, 0.09) &
  Gynaecological, Neuro-oncology and Orphan Cancer          & Student-t(0, 0.02, 1) & Normal(0, 0.10) &
  Gynaecology and Fertility                                 & Student-t(0, 0.02, 1) & Normal(0, 0.09) &
  Haematology                                               & Student-t(0, 0.05, 1) & Normal(0, 0.09) &
  Heart; Vascular                                           & Student-t(0, 0.08, 1) & Normal(0, 0.12) &
  Hepato-Biliary                                            & Student-t(0, 0.01, 1) & Normal(0, 0.09) &
  HIV/AIDS                                                  & Student-t(0, 0.02, 2) & Normal(0, 0.10) &
  Hypertension                                              & Student-t(0, 0.01, 2) & Normal(0, 0.06) &
  Incontinence                                              & Student-t(0, 0.06, 2) & Normal(0, 0.12) &
  Infectious Diseases                                       & Student-t(0, 0.02, 1) & Normal(0, 0.10) &
  Kidney and Transplant                                     & Student-t(0, 0.03, 1) & Normal(0, 0.09) &
  Lung Cancer                                               & Student-t(0, 0.06, 2) & Normal(0, 0.12) &
  Methodology                                               & Student-t(0, 0.10, 2) & Normal(0, 0.10) &
  Movement Disorders                                        & Student-t(0, 0.06, 2) & Normal(0, 0.07) &
  Multiple Sclerosis and Rare Diseases of the CNS           & Student-t(0, 0.07, 2) & Normal(0, 0.10) &
  Musculoskeletal                                           & Student-t(0, 0.01, 1) & Normal(0, 0.11) &
  Neonatal                                                  & Student-t(0, 0.02, 1) & Normal(0, 0.09) &
  Neuromuscular                                             & Student-t(0, 0.03, 1) & Normal(0, 0.08) &
  Oral Health                                               & Student-t(0, 0.03, 1) & Normal(0, 0.12) &
  Pain, Palliative and Supportive Care                      & Student-t(0, 0.11, 2) & Normal(0, 0.11) &
  Pregnancy and Childbirth                                  & Student-t(0, 0.01, 1) & Normal(0, 0.08) &
  Schizophrenia                                             & Student-t(0, 0.03, 1) & Normal(0, 0.10) &
  Sexually Transmitted Infections                           & Student-t(0, 0.04, 1) & Normal(0, 0.12) &
  Skin                                                      & Student-t(0, 0.03, 1) & Normal(0, 0.12) &
  Stroke                                                    & Student-t(0, 0.01, 1) & Normal(0, 0.06) &
  Tobacco Addiction                                         & Student-t(0, 0.04, 2) & Normal(0, 0.07) &
  Urology                                                   & Student-t(0, 0.04, 1) & Normal(0, 0.10) &
  Vascular                                                  & Student-t(0, 0.01, 1) & Normal(0, 0.12) &
  Work                                                      & Student-t(0, 0.06, 2) & Normal(0, 0.12) &
  Wounds                                                    & Student-t(0, 0.01, 1) & Normal(0, 0.09) &
  Cochrane                                                  & Student-t(0, 0.03, 1) & Normal(0, 0.10) "

  # update with manuscript
  #   Metabolic and Endocrine Disorders                         & Student-t(0, 0.00, 0) & Normal(0, 0.09) &
  #   Heart                                                     & Student-t(0, 0.00, 1) & Normal(0, 0.07) &
  #   Injuries                                                  & Student-t(0, 0.01, 0) & Normal(0, 0.14) &

  paper_logOR <- matrix(gsub("\n", "", strsplit(paper_logOR, "&")[[1]], fixed = TRUE), ncol = 3, byrow = T)
  paper_logRR <- matrix(gsub("\n", "", strsplit(paper_logRR, "&")[[1]], fixed = TRUE), ncol = 3, byrow = T)
  paper_RD    <- matrix(gsub("\n", "", strsplit(paper_RD,    "&")[[1]], fixed = TRUE), ncol = 3, byrow = T)

  for(i in 1:nrow(paper_logOR)){
    # test prior for the effect
    p1     <- prior_informed(paper_logOR[i,1], parameter = "effect", type = "logOR")
    p1pars <- as.numeric(strsplit(gsub(")", "", gsub("Student-t(", "", paper_logOR[i,2], fixed = TRUE), fixed = TRUE), ",")[[1]])
    expect_equal(p1$distribution, "t")
    expect_equal(p1$parameters, list("location" = p1pars[1], "scale" = p1pars[2], "df" = p1pars[3]))
    expect_equal(p1$truncation, list("lower" = -Inf, "upper" = Inf))

    # test prior for the heterogeneity
    if(gsub(" ", "", paper_logOR[i,3], fixed = TRUE) == "---"){
      expect_error(prior_informed(paper_logOR[i,1], parameter = "heterogeneity", type = "logOR"))
    }else{
      p2     <- prior_informed(paper_logOR[i,1], parameter = "heterogeneity", type = "logOR")
      p2pars <- as.numeric(strsplit(gsub(")", "", gsub("Inv-Gamma(", "", paper_logOR[i,3], fixed = TRUE), fixed = TRUE), ",")[[1]])
      expect_equal(p2$distribution, "invgamma")
      expect_equal(p2$parameters, list("shape" = p2pars[1], "scale" = p2pars[2]))
      expect_equal(p2$truncation, list("lower" = 0, "upper" = Inf))
    }
  }

  for(i in 1:nrow(paper_logRR)){
    # test prior for the effect
    p1     <- prior_informed(paper_logRR[i,1], parameter = "effect", type = "logRR")
    p1pars <- as.numeric(strsplit(gsub(")", "", gsub("Student-t(", "", paper_logRR[i,2], fixed = TRUE), fixed = TRUE), ",")[[1]])
    expect_equal(p1$distribution, "t")
    expect_equal(p1$parameters, list("location" = p1pars[1], "scale" = p1pars[2], "df" = p1pars[3]))
    expect_equal(p1$truncation, list("lower" = -Inf, "upper" = Inf))

    # test prior for the heterogeneity
    if(gsub(" ", "", paper_logRR[i,3], fixed = TRUE) == "---"){
      expect_error(prior_informed(paper_logRR[i,1], parameter = "heterogeneity", type = "logRR"))
    }else{
      p2     <- prior_informed(paper_logRR[i,1], parameter = "heterogeneity", type = "logRR")
      p2pars <- as.numeric(strsplit(gsub(")", "", gsub("Inv-Gamma(", "", paper_logRR[i,3], fixed = TRUE), fixed = TRUE), ",")[[1]])
      expect_equal(p2$distribution, "invgamma")
      expect_equal(p2$parameters, list("shape" = p2pars[1], "scale" = p2pars[2]))
      expect_equal(p2$truncation, list("lower" = 0, "upper" = Inf))
    }
  }

  for(i in 1:nrow(paper_RD)){
    # test prior for the effect
    p1     <- prior_informed(paper_RD[i,1], parameter = "effect", type = "RD")
    p1pars <- as.numeric(strsplit(gsub(")", "", gsub("Student-t(", "", paper_RD[i,2], fixed = TRUE), fixed = TRUE), ",")[[1]])
    expect_equal(p1$distribution, "t")
    expect_equal(p1$parameters, list("location" = p1pars[1], "scale" = p1pars[2], "df" = p1pars[3]))
    expect_equal(p1$truncation, list("lower" = -Inf, "upper" = Inf))

    # test prior for the heterogeneity
    if(gsub(" ", "", paper_RD[i,3], fixed = TRUE) == "---"){
      expect_error(prior_informed(paper_RD[i,1], parameter = "heterogeneity", type = "RD"))
    }else{
      p2     <- prior_informed(paper_RD[i,1], parameter = "heterogeneity", type = "RD")
      p2pars <- as.numeric(strsplit(gsub(")", "", gsub("Normal(", "", paper_RD[i,3], fixed = TRUE), fixed = TRUE), ",")[[1]])
      expect_equal(p2$distribution, "normal")
      expect_equal(p2$parameters, list("mean" = p2pars[1], "sd" = p2pars[2]))
      expect_equal(p2$truncation, list("lower" = 0, "upper" = Inf))
    }
  }

  expect_equal(print(prior_informed("cochrane", parameter = "effect", type = "logRR"), silent = TRUE),        "Student-t(0, 0.32, 3)")
  expect_equal(print(prior_informed("cochrane", parameter = "heterogeneity", type = "logRR"), silent = TRUE), "InvGamma(1.51, 0.23)")
})
