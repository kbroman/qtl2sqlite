context("Write and read probs")

test_that("write_probs and read_probs work", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    iron <- iron[,c(18,19,"X")]
    map <- insert_pseudomarkers(iron, step=0.5)
    pr <- calc_genoprob(iron, map, error_prob=0.002)

    file <- file.path(tempdir(), "iron.sqlite")

    write_probs(file, pr, map)

    pr18 <- read_probs(file, chr="18")
    expect_equal(pr18, pr[,18])

    prX <- read_probs(file, chr="X", pos=c(29.8, 32.8))
    expected <- pr[,"X"]
    expected$X <- expected$X[,,2:7]
    expect_equal(prX, expected)

    unlink(file)

})
