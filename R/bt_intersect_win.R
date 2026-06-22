bt_intersect_win <- function(region, record) {

  region <- data.table::as.data.table(region)
  record <- data.table::as.data.table(record)

  # 默认前三列为 chr/start/end
  names(region)[1:3] <- c("chr", "start", "end")
  names(record)[1:3] <- c("chr", "start", "end")

  # 保存原始列名
  region_col <- names(region)
  record_col <- names(record)

  # 添加id
  region$region_id <- seq_len(nrow(region))
  record$record_id <- seq_len(nrow(record))

  # 防止列名冲突
  names(region) <- paste0("region_", names(region))
  names(record) <- paste0("record_", names(record))

  # key
  data.table::setkey(
    region,
    region_chr,
    region_start,
    region_end
  )

  data.table::setkey(
    record,
    record_chr,
    record_start,
    record_end
  )

  # overlap
  hits <- data.table::foverlaps(
    record,
    region,
    by.x = c(
      "record_chr",
      "record_start",
      "record_end"
    ),
    by.y = c(
      "region_chr",
      "region_start",
      "region_end"
    ),
    nomatch = 0L
  )

  # 无 overlap
  if (nrow(hits) == 0) {
    return(data.frame())
  }

  # 真实交集区域
  hits$intersect_chr <- hits$region_chr

  hits$intersect_start <- pmax(
    hits$region_start,
    hits$record_start
  )

  hits$intersect_end <- pmin(
    hits$region_end,
    hits$record_end
  )

  hits$intersect_width <-
    hits$intersect_end -
    hits$intersect_start + 1

  return(as.data.frame(hits))
}
