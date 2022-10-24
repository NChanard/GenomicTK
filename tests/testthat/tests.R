GRange.gnr <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(c("chr1", "chr2"), c(3, 1)),
    ranges = IRanges::IRanges(c(1,201,251,1), end = c(200,250,330,100), names = letters[1:4]),
    strand = S4Vectors::Rle(BiocGenerics::strand(c("*")), 4),
    score = c(50,NA,100,30)
    )
GRange.gnr
chromSize.dtf = data.frame(c("chr1","chr2"),c(350,100))
binSize.num <- 100
binnedGRanges.gnr <- BinGRanges(
    gRange.gnr = GRange.gnr,
    chromSize.dtf=chromSize.dtf,
    binSize.num=binSize.num,
    method.chr ="mean",
    variablesName.chr_vec="score",
    na.rm=TRUE
)
binnedGRanges.gnr <- BinGRanges(
    gRange.gnr = GRange.gnr,
    chromSize.dtf=chromSize.dtf,
    binSize.num=binSize.num,
    method.chr ="mean",
    variablesName.chr_vec="score",
    na.rm=TRUE,
    cores.num=2
)

StrToGRanges("chr1:1-100:+")
StrToGRanges(c("chr1:1-100:+","chr2:400-500:-","chr1:10-50:*"))

GenomicSystem(1540,3)
GenomicSystem(1540,2)
GenomicSystem(10,2)
GenomicSystem(1000,2)
GenomicSystem(1000000,2)
GenomicSystem(1000000000,2)
GenomicSystem("1Gbp")
GenomicSystem("1Mbp")
GenomicSystem("1Kbp")
GenomicSystem("10Bp")

GRange_1.grn <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(c("chr1", "chr2", "chr1"), c(1, 3, 1)),
    ranges = IRanges::IRanges(101:105, end = 111:115, names = letters[1:5]),
    strand = S4Vectors::Rle(BiocGenerics::strand(c("-", "+", "*", "+")), c(1, 1, 2, 1)),
    score = 1:5
)
GRange_2.grn <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(c("chr1", "chr3"), c(1, 4)),
    ranges = IRanges::IRanges(106:110, end = 116:120, names = letters[6:10]),
    strand = S4Vectors::Rle(BiocGenerics::strand(c("*", "+", "-")), c(2, 1, 2)),
    score = 6:10
)
MergeGRanges(GRange_1.grn,GRange_2.grn)
GRange.lst = list(GRange_1.grn,GRange_2.grn)
MergeGRanges(GRange.lst)
MergeGRanges(GRange.lst, reduce.bln=TRUE, sort.bln=TRUE)

GRange.grn <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(c("chr1", "chr2", "chr1"), c(1, 3, 1)),
    ranges = IRanges::IRanges(101:105, end = 111:115, names = letters[1:5]),
    strand = S4Vectors::Rle(BiocGenerics::strand(c("-", "+", "*", "+")), c(1, 1, 2, 1)),
    seqinfo = c(chr1=200, chr2=300),
    score = 1:5
)
SeqEnds(GRange.grn)

StrToGRanges("chr1:1-100:+")
StrToGRanges(c("chr1:1-100:+","chr2:400-500:-","chr1:10-50:*"))
