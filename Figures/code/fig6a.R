library(grid)     ## Need to attach (and not just load) grid package
library(gtable)
library(pheatmap)

draw_colnames_50 <- function (coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 50, gp = gpar(...))
    return(res)}
# Modified pheatmap:::heatmap_motor
heatmap_motor <- function (matrix, border_color, cellwidth, cellheight, tree_col, 
    tree_row, treeheight_col, treeheight_row, filename, width, 
    height, breaks, color, legend, annotation_row, annotation_col, 
    annotation_colors, annotation_legend, annotation_names_row, 
    annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
    hjust_col, vjust_col, angle_col, fmat, fontsize_number, number_color, 
    gaps_col, gaps_row, labels_row, labels_col, ...) 
{
    lo = pheatmap:::lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), 
        ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, 
        treeheight_col = treeheight_col, treeheight_row = treeheight_row, 
        legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, 
        annotation_colors = annotation_colors, annotation_legend = annotation_legend, 
        annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, 
        main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
        fontsize_col = fontsize_col, angle_col = angle_col, gaps_row = gaps_row, 
        gaps_col = gaps_col, ...)
    res = lo$gt
    mindim = lo$mindim
    if (!is.na(filename)) {
        if (is.na(height)) {
            height = convertHeight(gtable_height(res), "inches", valueOnly = T)
        }
        if (is.na(width)) {
            width = convertWidth(gtable_width(res), "inches", valueOnly = T)
        }
        r = regexpr("\\.[a-zA-Z]*$", filename)
        if (r == -1) 
            stop("Improper filename")
        ending = substr(filename, r + 1, r + attr(r, "match.length"))
        f = switch(ending, pdf = function(x, ...) pdf(x, ...), 
            png = function(x, ...) png(x, units = "in", res = 300, 
                ...), jpeg = function(x, ...) jpeg(x, units = "in", 
                res = 300, ...), jpg = function(x, ...) jpeg(x, 
                units = "in", res = 300, ...), tiff = function(x, 
                ...) tiff(x, units = "in", res = 300, compression = "lzw", 
                ...), bmp = function(x, ...) bmp(x, units = "in", 
                res = 300, ...), stop("File type should be: pdf, png, bmp, jpg, tiff"))
        f(filename, height = height, width = width)
        gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, 
            border_color = border_color, tree_col = tree_col, 
            tree_row = tree_row, treeheight_col = treeheight_col, 
            treeheight_row = treeheight_row, breaks = breaks, 
            color = color, legend = legend, annotation_col = annotation_col, 
            annotation_row = annotation_row, annotation_colors = annotation_colors, 
            annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, 
            annotation_names_col = annotation_names_col, filename = NA, 
            main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
            fontsize_col = fontsize_col, hjust_col = hjust_col, 
            vjust_col = vjust_col, angle_col = angle_col, fmat = fmat, 
            fontsize_number = fontsize_number, number_color = number_color, 
            labels_row = labels_row, labels_col = labels_col, 
            gaps_col = gaps_col, gaps_row = gaps_row, ...)
        grid.draw(gt)
        dev.off()
        return(gt)
    }
    if (mindim < 3) 
        border_color = NA
    if (!is.na(main)) {
        elem = pheatmap:::draw_main(main, fontsize = 1.3 * fontsize, ...)
        res = gtable_add_grob(res, elem, t = 1, l = 3, name = "main", 
            clip = "off")
    }
    if (!pheatmap:::is.na2(tree_col) & treeheight_col != 0) {
        elem = pheatmap:::draw_dendrogram(tree_col, gaps_col, horizontal = T)
        res = gtable_add_grob(res, elem, t = 2, l = 3, name = "col_tree")
    }
    if (!pheatmap:::is.na2(tree_row) & treeheight_row != 0) {
        elem = pheatmap:::draw_dendrogram(tree_row, gaps_row, horizontal = F)
        res = gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
    }
    elem = pheatmap:::draw_matrix(matrix, border_color, gaps_row, gaps_col, 
        fmat, fontsize_number, number_color)
    res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
        name = "matrix")
    if (length(labels_col) != 0) {
        pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, 
            hjust_col = hjust_col, vjust_col = vjust_col, angle_col = angle_col, 
            ...)
        elem = do.call(pheatmap:::draw_colnames, pars)
        res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", 
            name = "col_names")
    }
    if (length(labels_row) != 0) {
        pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row, 
            ...)
        elem = do.call(pheatmap:::draw_rownames, pars)
        res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
            name = "row_names")
    }
    if (!pheatmap:::is.na2(annotation_col)) {
        converted_annotation = convert_annotations(annotation_col, 
            annotation_colors)
        elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
            gaps_col, fontsize, horizontal = T)
        res = gtable_add_grob(res, elem, t = 3, l = 3, clip = "off", 
            name = "col_annotation")
        if (annotation_names_col) {
            elem = pheatmap:::draw_annotation_names(annotation_col, fontsize, 
                horizontal = T)
            res = gtable_add_grob(res, elem, t = 3, l = 4, clip = "off", 
                name = "col_annotation_names")
        }
    }
    if (!pheatmap:::is.na2(annotation_row)) {
        converted_annotation = convert_annotations(annotation_row, 
            annotation_colors)
        elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
            gaps_row, fontsize, horizontal = F)
        res = gtable_add_grob(res, elem, t = 4, l = 2, clip = "off", 
            name = "row_annotation")
        if (annotation_names_row) {
            elem = pheatmap:::draw_annotation_names(annotation_row, fontsize, 
                horizontal = F, hjust_col = hjust_col, vjust_col = vjust_col, 
                angle_col = angle_col)
            res = gtable_add_grob(res, elem, t = 5, l = 2, clip = "off", 
                name = "row_annotation_names")
        }
    }
    annotation = c(annotation_col[length(annotation_col):1], 
        annotation_row[length(annotation_row):1])
    annotation = annotation[unlist(lapply(annotation, function(x) !pheatmap:::is.na2(x)))]
    if (length(annotation) > 0 & annotation_legend) {
        elem = pheatmap:::draw_annotation_legend(annotation, annotation_colors, 
            border_color, fontsize = fontsize, ...)
        t = ifelse(is.null(labels_row), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 6, b = 5, 
            clip = "off", name = "annotation_legend")
    }
    if (!pheatmap:::is.na2(legend)) {
        elem = pheatmap:::draw_legend(color, breaks, legend, fontsize = fontsize, 
            ...)
        t = ifelse(is.null(labels_row), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 5, b = 5, 
            clip = "off", name = "legend")
    }
    return(res)
}

# Modified pheatmap:::lo    
lo <- function (rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, 
    treeheight_col, treeheight_row, legend, annotation_row, annotation_col, 
    annotation_colors, annotation_legend, annotation_names_row, 
    annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
    angle_col, gaps_row, gaps_col, ...) 
{
    if (!is.null(coln[1]) | (!pheatmap:::is.na2(annotation_row) & annotation_names_row)) {
        if (!is.null(coln[1])) {
            t = coln
        }
        else {
            t = ""
        }
        tw = strwidth(t, units = "in", cex = fontsize_col/fontsize)
        if (annotation_names_row) {
            t = c(t, colnames(annotation_row))
            tw = c(tw, strwidth(colnames(annotation_row), units = "in"))
        }
        longest_coln = which.max(tw)
        gp = list(fontsize = ifelse(longest_coln <= length(coln), 
            fontsize_col, fontsize), ...)
        coln_height = unit(1, "grobheight", textGrob(t[longest_coln], 
            rot = angle_col, gp = do.call(gpar, gp))) + unit(10, 
            "bigpts")
    }
    else {
        coln_height = unit(5, "bigpts")
    }
    if (!is.null(rown[1])) {
        t = rown
        tw = strwidth(t, units = "in", cex = fontsize_row/fontsize)
        if (annotation_names_col) {
            t = c(t, colnames(annotation_col))
            tw = c(tw, strwidth(colnames(annotation_col), units = "in"))
        }
        longest_rown = which.max(tw)
        gp = list(fontsize = ifelse(longest_rown <= length(rown), 
            fontsize_row, fontsize), ...)
        rown_width = unit(1, "grobwidth", textGrob(t[longest_rown], 
            rot = 0, gp = do.call(gpar, gp))) + unit(10, "bigpts")
    }
    else {
        rown_width = unit(5, "bigpts")
    }
    gp = list(fontsize = fontsize, ...)
    if (!pheatmap:::is.na2(legend)) {
        longest_break = which.max(nchar(names(legend)))
        longest_break = unit(1.1, "grobwidth", 
            textGrob(as.character(names(legend))[longest_break], 
            gp = do.call(gpar, gp)))
        title_length = unit(1.1, "grobwidth", textGrob("Scale", 
            gp = gpar(fontface = "bold", ...)))
        legend_width = unit(12, "bigpts") + longest_break * 1.2
        legend_width = max(title_length, legend_width)
    }
    else {
        legend_width = unit(0, "bigpts")
    }
    if (is.na(main)) {
        main_height = unit(0, "npc")
    }
    else {
        main_height = unit(1.5, "grobheight", textGrob(main, 
            gp = gpar(fontsize = 1.3 * fontsize, ...)))
    }
    textheight = unit(fontsize, "bigpts")
    if (!pheatmap:::is.na2(annotation_col)) {
        annot_col_height = ncol(annotation_col) * (textheight + 
            unit(2, "bigpts")) + unit(2, "bigpts")
        t = c(as.vector(as.matrix(annotation_col)), colnames(annotation_col))
        annot_col_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
            gp = gpar(...))) + unit(12, "bigpts")
        if (!annotation_legend) {
            annot_col_legend_width = unit(0, "npc")
        }
    }
    else {
        annot_col_height = unit(0, "bigpts")
        annot_col_legend_width = unit(0, "bigpts")
    }
    if (!pheatmap:::is.na2(annotation_row)) {
        annot_row_width = ncol(annotation_row) * (textheight + 
            unit(2, "bigpts")) + unit(2, "bigpts")
        t = c(as.vector(as.matrix(annotation_row)), colnames(annotation_row))
        annot_row_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
            gp = gpar(...))) + unit(12, "bigpts")
        if (!annotation_legend) {
            annot_row_legend_width = unit(0, "npc")
        }
    }
    else {
        annot_row_width = unit(0, "bigpts")
        annot_row_legend_width = unit(0, "bigpts")
    }
    annot_legend_width = max(annot_row_legend_width, annot_col_legend_width)
    treeheight_col = unit(treeheight_col, "bigpts") + unit(5, 
        "bigpts")
    treeheight_row = unit(treeheight_row, "bigpts") + unit(5, 
        "bigpts")
    if (is.na(cellwidth)) {
        mat_width = unit(1, "npc") - rown_width - legend_width - 
            treeheight_row - annot_row_width - annot_legend_width
    }
    else {
        mat_width = unit(cellwidth * ncol, "bigpts") + length(gaps_col) * 
            unit(4, "bigpts")
    }
    if (is.na(cellheight)) {
        mat_height = unit(1, "npc") - main_height - coln_height - 
            treeheight_col - annot_col_height
    }
    else {
        mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * 
            unit(4, "bigpts")
    }
    gt = gtable(widths = unit.c(treeheight_row, rown_width,  
        mat_width, treeheight_row, legend_width, annot_legend_width), 
        heights = unit.c(main_height, treeheight_col, annot_col_height, 
            mat_height, coln_height), vp = viewport(gp = do.call(gpar, 
            gp)))
    cw = convertWidth(mat_width - (length(gaps_col) * unit(4, 
        "bigpts")), "bigpts", valueOnly = T)/ncol
    ch = convertHeight(mat_height - (length(gaps_row) * unit(4, 
        "bigpts")), "bigpts", valueOnly = T)/nrow
    mindim = min(cw, ch)
    res = list(gt = gt, mindim = mindim)
    return(res)
}

# Modified pheatmap:::draw_rownames      
draw_rownames <- function (rown, gaps, ...) 
{
    coord = pheatmap:::find_coordinates(length(rown), gaps)
    y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
    res = textGrob(rown, x = unit(-3, "bigpts"), y = y, vjust = 0.5, 
        hjust = 1, gp = gpar(...))
    return(res)
}

## 'Overwrite' default draw_rownames and draw_colnames with your own version 
assignInNamespace(x="draw_rownames", value=draw_rownames, ns="pheatmap")
assignInNamespace(x="lo", value=lo, ns="pheatmap")
assignInNamespace(x="heatmap_motor", value=heatmap_motor, ns="pheatmap")
assignInNamespace(x="draw_colnames", value="draw_colnames_50",
ns=asNamespace("pheatmap"))

discovery=read.table('data/discovery.set.counts.txt',header=T,comment.char = '!')
abbs=unique(read.table('data/GWAS.abbv.txt',header=T)[,1:2])
rownames(abbs)=abbs[,1]

dat.signif=read.table('data/GWAScoloc.signif_fdr001.tissuespecific.txt',header=T,row.names=1,sep='\t',colClasses = 'character')
dat=read.table('data/GWAScoloc.scaled.props.tissuespecific.txt',header=T,row.names=1)
rownames(dat)=discovery$Tissue
colnames(dat)=as.character(abbs[colnames(dat),2])
breaksList = seq(round(-max(abs(dat)),1),round(max(abs(dat)),1), by = 0.05)
pheatmap(dat,color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),fontsize_col=13,fontsize_row = 13,fontsize = 13,display_numbers = dat.signif,cluster_rows = F,angle_col=90,hjust_col=1,vjust_col=1,breaks=breaksList,labels_col=as.character(colnames(dat)))
#600x500
