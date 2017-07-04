N.pair <- readRDS("~/clb/ryckman/wes/data/N-pair.rds")
N.trio <- readRDS("~/clb/ryckman/wes/data/N-trio.rds")
vData <- readRDS("~/clb/ryckman/wes/data/vData.rds")

geneTest("SAMD11", N.pair, N.trio, vData)
geneTest("SAMD11", N.pair)

geneSummary("SAMD11", N.pair, N.trio)

rvis("SAMD11")
rvis(c("SAMD11", "ISG15"))
