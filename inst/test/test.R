library(epcot)
N2 <- epcot_example$N2
N3 <- epcot_example$N3
vData <- epcot_example$vData
res <- bayesTest(N2, N3, vData)

expect_equal(nrow(N2), 1000)
expect_equal(nrow(N3), 1000)
expect_equal(nrow(vData), 1000)

geneTest("SAMD11", N2, N3, vData=vData)
geneTest("SAMD11", N2, vData=vData)

geneSummary("SAMD11", N2, N3, vData=vData)

rvis("SAMD11")
rvis(c("SAMD11", "ISG15"))
