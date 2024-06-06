# https://blog.r-hub.io/2019/12/12/internal-functions/

#' Check arguments
#'
#' Check if arguments are properly formatted and fit the ldf
#' @noRd
check.arguments = function(ldf = parent.frame()$ldf,
                               value.col = parent.frame()$value.col,
                               group.by = parent.frame()$group.by,
                               fill.by = parent.frame()$fill.by,
                               copy.by = parent.frame()$copy.by,
                               exclude.by = parent.frame()$exclude.by,
                               apply.to = parent.frame()$apply.to,
                               apply.not.to = parent.frame()$apply.not.to,
                               test = parent.frame()$test,
                               test.not = parent.frame()$test.not,
                               excl = parent.frame()$excl,
                               excl.not = parent.frame()$excl.not,
                               view = parent.frame()$view,
                               view.not = parent.frame()$view.not,
                               off.limit.colnames = NULL){
  if(!is.data.frame(ldf)) stop("'ldf' must be a dataframe")
  cn.df = colnames(ldf)
  if(ncol(ldf) < 2) return("'ldf' must have one column containing the values and at least one identifier-column used to group them.")
  if(length(cn.df) != length(unique(cn.df))) return("Colnames of 'ldf' must be unique")

  if(any(off.limit.colnames %in% cn.df)) return(paste0("[",paste(off.limit.colnames[off.limit.colnames %in% cn.df], collapse = ","), "] may not be used as column name"))
  if(!value.col %in% cn.df) return(paste0("[",value.col,"] required by 'value.col' but not present in the colnames"))

  groupers = list(group.by = group.by, fill.by = fill.by, copy.by = copy.by, exclude.by = exclude.by)
  for(i in names(groupers)){
    if(value.col %in% groupers[[i]]) return(paste("The value-col cannot be defined in '",i,"'"))
    if(any(!groupers[[i]] %in% cn.df)) return(paste0("[",paste(groupers[[i]][!groupers[[i]] %in% cn.df], collapse = ","), "] required by '",i,"' but not present in the colnames of 'ldf'"))
    if(sum(is.na(ldf[,groupers[[i]]])) > 0) return(paste0("The identifier-columns defined by '",i,"' can not contain NA"))
    if(i == 'exclude.by' & !is.null(exclude.by)) if(any(!exclude.by %in% group.by)) return("'exclude.by' can only use columns from 'group.by' (grouping used for exclusion cannot be more specific than grouping used for testing)")
    if(any(copy.by %in% group.by)) return("'copy.by' and 'group.by' must use different columns")
  }

  selectors = list(apply.to = apply.to, apply.not.to = apply.not.to, test = test, test.not = test.not, excl = excl, excl = excl.not, view = view, view.not = view.not)
  for(i in names(selectors)){
    if(!is.null(selectors[[i]])){
      if(!is.list(selectors[[i]])) return(paste0("'",i,"' must be a list"))
      if(any(is.null(names(selectors[[i]])))) return(paste0("'",i,"' must have element names corresponding to column names of 'ldf'"))
      if(any(!names(selectors[[i]]) %in% colnames(ldf))) return(paste0("[",paste(names(selectors[[i]])[!names(selectors[[i]]) %in% colnames(ldf)], collapse = ','), "] required by '",i,"', but not present in the colnames of 'ldf'"))
      for(j in names(selectors[[i]])){
        if(any(!selectors[[i]][[j]] %in% ldf[[j]])) return(paste0(paste(selectors[[i]][[j]][!selectors[[i]][[j]] %in% ldf[[j]]], collapse = ','), " required by '",i,"', but not present in column [",j,"] in 'ldf'"))
        if(i %in% c('apply.not.to','test.not','excl.not','view.not'))
          if(all(ldf[[j]] %in% selectors[[i]])) return(paste0("All identifiers of column [",i,"] are present in '",i,"' and no rows would be left to be used."))}}
  }
}



#' Select rows from an ldf
#'
#' Use apply.to and apply.not.to to select or deselect rows
#' @param ldf A dataframe
#' @param apply.to,apply.not.to Named lists of vectors, with the name corresponding to the colnames of the ldf and the vectors containing the identifiers used to (un)select rows from ldf that contain the identifier(s) in the specified columns.
#' @noRd
apply.to.not.to = function(ldf, apply.to, apply.not.to, return.bool = F){
  bool = rep(T,nrow(ldf))
  if(!is.null(apply.to))     for (i in names(apply.to))     if(i %in% colnames(ldf)) bool[!ldf[[i]] %in% apply.to[[i]]] = F
  if(!is.null(apply.not.to)) for (i in names(apply.not.to)) if(i %in% colnames(ldf)) bool[ ldf[[i]] %in% apply.not.to[[i]]] = F
  # if(!is.null(apply.to))     for (i in names(apply.to))     if(i %in% colnames(ldf)) ldf = ldf[ ldf[[i]] %in% apply.to[[i]],]
  # if(!is.null(apply.not.to)) for (i in names(apply.not.to)) if(i %in% colnames(ldf)) ldf = ldf[!ldf[[i]] %in% apply.not.to[[i]],]
  if(return.bool) return(bool) else return(ldf[bool,])}


#' A function to create example trial-level data
#'
#' A function to create example trial-level local-global data from some default experiment to be used for example calculations
#' @param mean.RT,sd.RT Scalar, overall mean and standard deviation values being passed to rnorm() to create example reactiont times. Any negative values generated by rnorm() will be made positive. Default to 400 and 150, respectively.
#' @param mean.accuracy Scalar, overall mean accuracy used to create example accuracies. Defaults to 0.9.
#' @param condition,shuffle A vector of trial conditions, shuffled (T) or not shuffled (F) for each participant and each session. Is being used to determine the number of trials for each participant, repeated for each session. Defaults to c(rep('A',5), rep('B',5)).
#' @param id A vector of participant identifiers. Defaults to '1:10'
#' @param session A vector of session identifiers. Defaults to '1:2'
#' @keywords example, data
#' @export
#' @examples example.trial.data()
example.trial.data = function(seed = NULL, mean.RT = 400, sd.RT = 150, mean.accuracy = 0.9, condition = c(rep('A',5), rep('B',5)), shuffle = T, id = 1:8, session = 1:2){
  set.seed(seed = seed)
  # Create a long dataframe of example trial data for 10 participants in a 2-session experiment with 5 A-trials and 5 B-trials trials, half of them congruent and half incongruent
  tr = expand.grid(list(measure = c('RT','ACC'), condition = condition, id = id, session = session))
  tr$trial = rep(1:length(condition), each = 2)
  tr$value = NA
  tr$value[seq(2,nrow(tr),2)] = sample(c(rep(0,round(1-mean.accuracy, 3) * 1000), rep(1,round(mean.accuracy, 3) * 1000)), nrow(tr)/2, replace = T)
  tr$value[seq(1,nrow(tr),2)] = abs(rnorm(nrow(tr)/2, mean = mean.RT, sd = sd.RT))
  if(shuffle) tr$condition = rep(as.vector(apply(matrix(tr[tr$measure == 'RT','condition'], nrow = length(condition)), 2, sample)), each = 2) # If required, shuffle the conditions for each participant and session
  return(tr[,c('session','id','trial','condition','measure','value')])
}




#' Fill up missing values in an ldf
#'
#' Fill missing rows to make sure every combination of the grouping-identifiers, measure-identifiers, and observations, exists (measure and observation can be disabled)
#' @param ldf A long format dataframe with one column holding the values and at least one column holding identifiers that can be used to group the values. Compare example.trial.data().
#' @param value.col Name of the column holding the values in the ldf. Defaults to 'value'.
#' @param fill.by Identifier-columns that are used to create groupings, which are filled up. Defaults to NULL, which means all columns besides 'value.col' are used.
#' @param notify Provides progress information. Defaults to TRUE.
#' @importFrom dplyr '%>%'
#' @keywords data
#' @export
#' @examples
#' # Introduce some missing rows into the example data
#' ldf = example.trial.data(); ldf = ldf[sort(sample(1:nrow(ldf), 0.6*nrow(ldf))),]
#' # Fill the rows of all missing groupings
#' ldf2 = fill.missing.values(ldf = ldf)
#' # Note that not specifying 'fill.by', as done before, leads to non-desired output in the example dataset, because each value is already fully specified with the columns session,id,trial,measure. The column condition introduces over-specification, since when trial 1 for participant 1 in session 1 has EITHER condition A OR condition B. If this row 1,1,1 is missing, it needs to be added but with condition NA, because the condition of the trial is not known. Therefore, condition should not be included in 'fill.by'.
#' ldf2 = fill.missing.values(ldf = ldf, fill.by = c('id','session','trial','measure'))
fill.missing.values = function(ldf, fill.by = NULL, notify = T, value.col = 'value'){

  # Check if args match the ldf
  stp = check.arguments(); if(!is.null(stp)) stop(stp)

  # If fill.by is not provided, use all identifier-columns
  if(is.null(fill.by)) fill.by = colnames(ldf)[!colnames(ldf) %in% value.col]

  # Get the original number of rows and size of groupings
  nrow.ldf = nrow(ldf)
  group.size = ldf %>% dplyr::group_by(dplyr::across(dplyr::all_of(fill.by))) %>% dplyr::group_size()

  # Create a ldf that contains all possible combinations of the identifiers in the fill.by cols
  ldf2 = expand.grid(lapply(ldf[,fill.by], unique))

  # Merge ldf and ldf2, rows of ldf2 and cols of ldf are maintained, nonexistant values become NA
  ldf = suppressMessages({dplyr::left_join(ldf2, ldf)})

  # Notifications and stats
  if(notify & ncol(ldf) > length(fill.by) + 1) message("Column [", paste(colnames(ldf)[!colnames(ldf) %in% c(fill.by,value.col)], collapse = ','),"] not specified in 'fill.by'. NAs are introduced in new rows.")
  if(notify) message(nrow(ldf) - nrow.ldf, " rows added to ",nrow.ldf," existing rows.")
  if(length(unique(ldf %>% dplyr::group_by(dplyr::across(dplyr::all_of(fill.by))) %>% dplyr::group_size())) > 1) warning("After filling, group sizes are still not identical. Possibly you missed to specify a relevant identifier column.")

  return(ldf)
}



#' Exclude groupings of an ldf, if they fulfil a certain criterion
#'
#' Groups an ldf by certain columns, calculates a freely definable statistic for each grouping, compares the result against a certain threshold, and trims all values in each grouping that crosses the threshold by overwriting them with NA
#' @param ldf A long format dataframe with one column (assumed name: 'value') holding the values and at least one column holding identifiers that can be used to group the values. See example.trial.data() for a correctly formatted long dataframe.
#' @param value.col Name of the column holding the values in the ldf. Defaults to 'value'.
#' @param group.by Identifier-columns that are used to create groupings, to which the test is applied.
#' @param exclude.by Identifier-columns that are used to group the data for exclusion. Can not contain names that are not in 'group.by' (i.e. can not be more specific). Defaults to NULL, which means names will be copied from group.by
#' @param test,test.not,excl,excl.not A named list of vectors to (un)select rows in the testing/exclusion phase. The element names correspond to column names and the vectors contain the identifiers which should (not) be considered when testing/excluding groupings. Default to NULL, which means no (un)selection is applied to any column.
#' @param statistic Name of a function used to calculate a summary statistic for the values of each grouping. Must take a numeric vector and output a numeric scalar. Defaults to 'mean'. 'na.rm = T' is passed as additional argument.
#' @param criterion Function that defines when a grouping statistic (output of the 'statistic') should be excluded. Defaults to 'function(S) return(S <= .5)', which means all groupings with a summary-statistic lower equal 0.5 ought to be excluded. The function must accept a numeric scalar/vector as input and return a logical scalar/vector.
#' @param report.by Identifier-columns that are used to create groupings, for which statistics are printed when notify is TRUE. Defaults to NULL, which means one value for the entire ldf is returned.
#' @param return.report,return.stats Adds the trimming report/summary statistics to the ldf in a list. Defaults to FALSE.
#' @param notify Provides progress information. Defaults to TRUE.
#' @importFrom dplyr '%>%'
#' @keywords exclusion
#' @export
#' @examples
#' # Exclude each participant in each session who has a mean ACC smaller or equal 0.8. 'exclude.by' is not defined, therefore a grouping's test is only applied to itself. 'excl' is not defined, therefore RT values are excluded as well, not only the ACC-values which were tested
#' ldf = exclude.groupings(ldf = example.trial.data(seed = 1), group.by = c('id','session'), test = list(measure = 'ACC'), criterion = function(S) return(S <= 0.8))
#' # Do the same, but now make the exclusion-grouping less fine-grained than the test-grouping: By setting exclude.by = 'id', a participant's values are excluded in both sessions, if the participant fulfils the exclusion criteria in either or both sessions (a test of the 'id*session'-sub-grouping is projected to the exclusion-super-grouping 'id')
#' ldf = exclude.groupings(ldf = example.trial.data(), group.by = c('id','session'), exclude.by = 'id', test = list(measure = 'ACC'), criterion = function(S) return(S <= 0.8))
#' # Do the same, but now exempt session 2 from exclusion. Note that it is not exempt from testing, which means a participant who fulfils the criterion in session 2 is excluded in session 1, but session 2 is exempt
#' ldf = exclude.groupings(ldf = example.trial.data(), group.by = c('id','session'), exclude.by = 'id', test = list(measure = 'ACC'), excl.not = c(session = '2'), criterion = function(S) return(S <= 0.8))
#' # Back to the basic exclusion, but while testing ACC, only apply exclusion to RT. 'report.by' won't show
#' ldf = exclude.groupings(ldf = example.trial.data(), group.by = c('id','session'), test = list(measure = 'ACC'), excl = list(measure = 'RT'), criterion = function(S) return(S <= 0.8), report.by = c('measure','session'))
#' stat = trim.statistics(ldf = ldf, report.by = c('measure','session'), print.stats = T)
#' # Exclude participants entire ACC data, if their median RT in condition A was higher than 1000 or lower than 400
#' ldf = exclude.groupings(ldf = example.trial.data(), group.by = c('id','session','condition'), exclude.by = 'id', test = list(measure = 'RT', condition = 'A'), excl = list(measure = 'ACC'), statistic = 'median', criterion = function(S) return(S < 400 | S > 1000), report.by = 'measure')
exclude.groupings = function(ldf, group.by, exclude.by = NULL, test = NULL, test.not = NULL, excl = NULL, excl.not = NULL, statistic = 'mean', criterion = function(S) return(S <= 0.5), value.col = 'value', report.by = NULL, return.report = F, return.stats = F, notify = T){

  # Check whether the dataframe and the args have the right prorperties and fit together
  stp = check.arguments(off.limit.colnames = c('grouping.statistic','exclusion.column')); if(!is.null(stp)) stop(stp)
  # Check if the function defined for criterion only takes one argument
  if(length(formals(criterion)) > 1) stop("THe function defined as 'criterion' can only take a single argument (the statistic calculated by 'statistic')")

  # If exclude.by is NULL, it is the same as group.by
  if(is.null(exclude.by)) exclude.by = group.by

  # Group the ldf and create an original version (ldf0) and a version to work with (ldf2)
  ldf0 = ldf = ldf2 = ldf %>% dplyr::group_by(dplyr::across(dplyr::all_of(group.by)))

  # Remove the rows which should not affect the exclusion from the working-data ldf2
  ldf2 = apply.to.not.to(ldf = ldf2, apply.to = test, apply.not.to = test.not)

  # Apply the 'statistic' to the values of the groupings
  stat = suppressMessages({ldf2 %>% dplyr::summarise(grouping.statistic = do.call(match.fun(statistic), list(!!dplyr::sym(value.col), na.rm = T)), .groups = 'keep')})
  # Find groupings whose stat meets the criterion
  stat$exclusion.column = criterion(stat$grouping.statistic)
  # Reduce the exclusion template with 'exclude.by'. If exclude.by has a lower 'resolution' it merges groupings that were separated by 'group.by'. If any of these merged groupings were TRUE, the entire exclusion-grouping becomes TRUE
  stat2 = suppressMessages({stat %>% dplyr::group_by(across(all_of(exclude.by))) %>% dplyr::summarise(exclusion.column = any(exclusion.column), .groups = 'keep')})

  # Remove the rows which should not be affected by exclusion from the return-data ldf and the exclusion template stat2
  ldf = apply.to.not.to(ldf = ldf, apply.to = excl, apply.not.to = excl.not)
  stat2 = apply.to.not.to(ldf = stat2, apply.to = excl, apply.not.to = excl.not)

  # TODO is the second statement really correct? It seemed once that it just reiterates stat
  if(notify){
    message(sum(stat$exclusion.column)," of ",nrow(stat)," tested (and ",dplyr::n_groups(ldf0)," total) groupings (by [",paste(group.by,collapse = ','),"]) fulfil the exclusion criterion")
    message(sum(stat2$exclusion.column)," of ",dplyr::n_groups(ldf %>% dplyr::group_by(across(all_of(exclude.by))))," selected (and ",dplyr::n_groups(ldf0 %>% dplyr::group_by(across(all_of(exclude.by))))," total) groupings (by [",paste(exclude.by,collapse = ','),"]) are excluded")}

  # Join the exclusion template stat2 with the ldf. If rows in stat3 do not have an equivalent in ldf, they are discarded (in doing this, relevant exclusion-information is not lost, because in the reduce-step from stat to stat2, all relevant information from test-(sub-)groupings was extracted and fed into exclusion-(super-)groupings. Information that is dropped now, is per user-requirement to be neglected)
  ldf = suppressMessages({dplyr::left_join(ldf, stat2)})
  ldf$exclusion.column[is.na(ldf$exclusion.column)] = F
  # Join the ldf with ldf0
  ldf = suppressMessages({dplyr::left_join(ldf0, ldf)})
  # NA those values where excl is TRUE
  ldf[[value.col]][dplyr::coalesce(ldf$exclusion.column, F)] = NA
  ldf = ldf %>% dplyr::select(-exclusion.column)

  # Notification of trimming
  report = trim.statistics(ldf = ldf, ldf0 = ldf0, report.by = report.by, apply.to = excl, apply.not.to = excl.not, print.stats = notify, value.col = value.col)

  # Nicer names for stats
  colnames(stat)[colnames(stat) %in% 'grouping.statistic'] = statistic
  colnames(stat)[colnames(stat) %in% 'exclusion.column'] = 'excl'

  # Return with or without report and stats
  rt = list(ldf = as.data.frame(ldf))
  if(return.report) rt = c(rt, report = report)
  if(return.stats) rt = c(rt, stats = stats)
  if(return.report | return.stats) return(rt) else return(as.data.frame(ldf))
}


#' Exclude groupings of an ldf, if they fulfil a certain criterion
#'
#' Groups an ldf by certain columns, calculates a freely definable statistic for each grouping, compares the result against a certain threshold, and trims all values in each grouping that crosses the threshold by overwriting them with NA
#' @param ldf A long format dataframe with one column (assumed name: 'value') holding the values and at least one column holding identifiers that can be used to group the values. See example.trial.data() for a correctly formatted long dataframe.
#' @param value.col Name of the column holding the values in the ldf. Defaults to 'value'.
#' @param group.by Identifier-columns that are used to create groupings, to which the test is applied.
#' @param exclude.by Identifier-columns that are used to group the data for exclusion. Can not contain names that are not in 'group.by' (i.e. can not be more specific). Defaults to NULL, which means names will be copied from group.by
#' @param test,test.not,excl,excl.not A named list of vectors to (un)select rows in the testing/exclusion phase. The element names correspond to column names and the vectors contain the identifiers which should (not) be considered when testing/excluding groupings. Default to NULL, which means no (un)selection is applied to any column.
#' @param statistic Name of a function used to calculate a summary statistic for the values of each grouping. Must take a numeric vector and output a numeric scalar. Defaults to 'mean'. 'na.rm = T' is passed as additional argument.
#' @param thresholds Numeric vector of length two, forming the lower and upper bound. Groupings whose summary statistic falls on or outside these bounds are excluded. NA are interpreted as Inf. Defaults to c(Inf,Inf).
#' @param flip Flips the logic of the threshold, and excludes all values that fall on or inside the upper and lower bound. Defauts to FALSE.
#' @param incl.thresholds Makes the threshold-values inclusive. Defaults to TRUE.
#' @param report.by Identifier-columns that are used to create groupings, for which statistics are printed when notify is TRUE. Defaults to NULL, which means one value for the entire ldf is returned.
#' @param return.report,return.stats Adds the trimming report/summary statistics to the ldf in a list. Defaults to FALSE.
#' @param notify Provides progress information. Defaults to TRUE.
#' @importFrom dplyr '%>%'
#' @keywords exclusion
#' @export
#' @examples
#' # Exclude each participant in each session who has a mean ACC smaller or equal 0.8. 'exclude.by' is not defined, therefore a grouping's test is only applied to itself. 'excl' is not defined, therefore RT values are excluded as well, not only the ACC-values which were tested
#' ldf = exclude.groupings(ldf = example.trial.data(seed = 1), group.by = c('id','session'), test = list(measure = 'ACC'), criterion = function(S) return(S <= 0.8))
#' # Do the same, but now make the exclusion-grouping less fine-grained than the test-grouping: By setting exclude.by = 'id', a participant's values are excluded in both sessions, if the participant fulfils the exclusion criteria in either or both sessions (a test of the 'id*session'-sub-grouping is projected to the exclusion-super-grouping 'id')
#' ldf = exclude.groupings(ldf = example.trial.data(), group.by = c('id','session'), exclude.by = 'id', test = list(measure = 'ACC'), criterion = function(S) return(S <= 0.8))
#' # Do the same, but now exempt session 2 from exclusion. Note that it is not exempt from testing, which means a participant who fulfils the criterion in session 2 is excluded in session 1, but session 2 is exempt
#' ldf = exclude.groupings(ldf = example.trial.data(), group.by = c('id','session'), exclude.by = 'id', test = list(measure = 'ACC'), excl.not = c(session = '2'), criterion = function(S) return(S <= 0.8))
#' # Back to the basic exclusion, but while testing ACC, only apply exclusion to RT. 'report.by' won't show
#' ldf = exclude.groupings(ldf = example.trial.data(), group.by = c('id','session'), test = list(measure = 'ACC'), excl = list(measure = 'RT'), criterion = function(S) return(S <= 0.8), report.by = c('measure','session'))
#' stat = trim.statistics(ldf = ldf, report.by = c('measure','session'), print.stats = T)
#' # Exclude participants entire ACC data, if their median RT in condition A was higher than 1000 or lower than 400
#' ldf = exclude.groupings(ldf = example.trial.data(), group.by = c('id','session','condition'), exclude.by = 'id', test = list(measure = 'RT', condition = 'A'), excl = list(measure = 'ACC'), statistic = 'median', criterion = function(S) return(S < 400 | S > 1000), report.by = 'measure')
exclude.groupings.by_interval = function(ldf, group.by, exclude.by = NULL, test = NULL, test.not = NULL, excl = NULL, excl.not = NULL, statistic = 'mean', thresholds = c(Inf, Inf), flip = F, incl.thresholds = T, value.col = 'value', report.by = NULL, return.report = F, return.stats = F, notify = T){

  # Check the thresholds argument
  if(length(thresholds) != 2) stop("'thresholds' must be of length 2")
  thresholds[is.na(thresholds)] = Inf

  # Create the criterion-function, based on the specification of incl.thresholds and flip
  if(incl.thresholds){
    if(!flip){
      # Exclude everything outside the defined thresholds
      interval.criterion = function(S) return(S <= thresholds[1] | S >= thresholds[2])
      interval.criterion = function(S) return(S >= thresholds[1] & S <= thresholds[2])}
    } else {
      # Exclude everything inside the defined thresholds
      if(!flip){
        interval.criterion = function(S) return(S < thresholds[1] | S > thresholds[2])
        interval.criterion = function(S) return(S > thresholds[1] & S < thresholds[2])}
    }
  # Pass all arguments to the exclude.groupings base function and return its output
  return(excludr::exclude.groupings(ldf = ldf, group.by = group.by, exclude.by = exclude.by, test = test, test.not = test.not, excl = excl.not, statistic = statistic, criterion = interval.criterion, value.col = value.col, report.by = report.by, return.report = return.report, return.stats = return.stats, notify = notify))}



trim.values.by_interval = function(ldf, excl = NULL, excl.not = NULL, threshold = c(Inf,Inf), flip = F, incl.thresholds, value.col = 'value'){
  if(length(thresholds) > 2) stop("'thresholds' can not be longer than 2")
  if(length(thresholds) == 1) thresholds = c(thresholds,thresholds)
  thresholds[is.na(thresholds)] = Inf

  # Create a bool where the rows that may be excluded are T
  bool1 = excludr::apply.to.not.to(ldf = ldf, apply.to = excl, apply.not.to = excl.not, return.bool = T)
  # Create a bool for which values fulfill the threshold criteria
  if(incl.thresholds){
    if(!flip){
      bool2 = ldf[[value.col]] <= threshold[1] | ldf[[value.col]] >= threshold[2]
      bool2 = ldf[[value.col]] >= threshold[1] | ldf[[value.col]] <= threshold[2]}
  } else {
    if(!flip){
      bool2 = ldf[[value.col]] < threshold[1] & ldf[[value.col]] > threshold[2]
      bool2 = ldf[[value.col]] > threshold[1] & ldf[[value.col]] < threshold[2]}
  }
  # Trim the values with the bools
  ldf[[value.col]][bool1 & bool2] = NA
  return(ldf)
}



#' Trim the lowest and/or highest values in groupings of an ldf
#'
#' Groups an ldf by certain columns and finds the x highest and or lowest values in each grouping, to overwrite them with NA. Tied values can be trimmed together or randomly.
#' @param ldf A long format dataframe with one column holding the values and at least one column holding identifiers that can be used to group the values. Compare example.trial.data().
#' @param value.col Name of the column holding the values in the ldf. Defaults to 'value'.
#' @param group.by Identifier-columns that are used to create groupings, to which the test is applied.
#' @param apply.to,apply.not.to A named list of vectors to (un)select rows. The names of the columns are specified by the element names, and the identifiers in this column to which the trimming should (not) be applied are specified by the vector. Default to NULL, which means no (un)selection is applied to any column.
#' @param report.by Identifier-columns that are used to create groupings, for which statistics are returned.
#' @param trim Proportion or number of values in a grouping that should be trimmed from bottom and top. Values < 1 are interpreted as proportion of observations, values >= 1 as number of observations. If only one value is provided, it will be used for bottom and top. Defaults to c(0.025,0.025), i.e. 95&#37; interval
#' @param trim.ties Tied values around the edge of the trim-cutoff defined by 'trim' are trimmed together, even when doing so increases the share of trimmed values beyond the specification of trim. If FALSE, tied values to be trimmed are sampled, which can lead to unexpected behaviour when the ties around the bot-cutoff and the top-cutoff are the same (e.g. when a grouping has only a single unique value). Defaults to TRUE.
#' @param warn.ties Issues a warning when there are ties in any grouping defined by 'group.by'. Defaults to TRUE.
#' @param view.ties Returns a dataframe with descriptions of the ties instead of the trimmed ldf. Defaults to FALSE.
#' @param ignore.na Ignores NA-values when calculating the number of values to be trimmed from the trim-proportion and the size of the grouping. Defaults to TRUE.
#' @param report.by Identifier-columns that are used to create groupings, for which statistics are printed when notify is TRUE. Defaults to NULL, which means one value for the entire ldf is returned.
#' @param return.report Adds the trimming report to the trimmed ldf in a list. Defaults to FALSE.
#' @param notify Provides progress information. Defaults to TRUE.
#' @importFrom dplyr '%>%'
#' @keywords trimming
#' @export
#' @examples
#' # Trim the bottom and top 0.1 RT values for each participant in each session separately
#' ldf = trim.by_interval(ldf = example.trial.data(), group.by = c('id', 'session', 'measure'), apply.to = list(measure = 'RT'), trim = 0.1)
#' # Do the same again using the already trimmed ldf. Note that no values were trimmed and a warning was issued: ldf now only has 8 non-NA values per grouping and since NAs are ignored per default, trim = 0.1 results in a requirement to trim 0.8 values, which is rounded down to 0.
#' ldf2 = trim.by_interval(ldf = ldf, group.by = c('id', 'session'), apply.to = list(measure = 'RT'), trim = 0.1)
#' # When ignore.na = FALSE or the proportion to be trimmed is higher, this does not happen
#' ldf2 = trim.by_interval(ldf = ldf, group.by = c('id', 'session'), apply.to = list(measure = 'RT'), trim = 0.1, ignore.na = F)
#' ldf2 = trim.by_interval(ldf = ldf, group.by = c('id', 'session'), apply.to = list(measure = 'RT'), trim = 0.13)
#' # Trim the four highest RT values, when additionally grouping by 'condition' and excluding participant 1 and 2.
#' ldf = trim.by_interval(ldf = example.trial.data(), group.by = c('id', 'session', 'measure'), apply.to = list(measure = 'RT'), apply.not.to = list(id = 1:2), trim = c(0,4))
#' # Trim the highest and lowest observation and report the trimming for each measure: All ACC-values are trimmed, because it only contains the values 0 and 1, and since trim.ties = T is default, all 0s and 1s are trimmed when 0 and 1 are the highest and lowest value.
#' ldf = trim.by_interval(ldf = example.trial.data(), group.by = c('id', 'session', 'measure'), report.by = 'measure', trim = 1)
#' # Set 'trim.ties' to FALSE to avoid this behavior and trim strictly the amount of values defined by 'trim'. Values to be trimmed are sampled randomly from the pool of relevant tied values. If the tied value affecting the upper and lower threshold is the same (as can be the case here for ACC, with all 1s), it can happen that the mechanism for the upper and lower threshold sample the same row and thus fewer values than required by 'trim' are trimmed.
#' ldf = trim.by_interval(ldf = example.trial.data(), group.by = c('condition','id', 'session','measure'), report.by = 'measure', trim = 1, trim.ties = F)
trim.values.by_quantile = function(ldf, group.by, apply.to = NULL, apply.not.to = NULL, trim = 0.025, trim.ties = T, warn.ties = T, view.ties = F, ignore.na = T, report.by = NULL, return.report = F, notify = T, value.col = 'value'){

  # Check if args match the ldf
  stp = check.arguments(); if(!is.null(stp)) stop(stp)
  # Check if all groupings have the same size
  if(length(unique(dplyr::group_size(ldf %>% dplyr::group_by(dplyr::across(dplyr::all_of(group.by)))))) > 1) warning("Not all groups as defined by group.by have the same size. Unique grouping sizes are: ",paste(unique(dplyr::group_size(ldf %>% dplyr::group_by(dplyr::across(dplyr::all_of(group.by))))), collapse = ', '))

  # Ensure 'trim' has length 2
  if(length(trim) > 2) stop("'trim' can contain 1 or two 2 values")
  if(length(trim) == 1) trim = rep(trim,2)

  # Group the ldf and create a second ldf, in which the trimming is done
  ldf = ldf %>% dplyr::group_by(dplyr::across(dplyr::all_of(group.by)))
  dupl = ldf %>% select(-!!dplyr::sym(value.col)) %>% dplyr::group_by_all() %>% dplyr::group_size()
  if(any(dupl > 1)) stop(sum(dupl-1)," of ",length(dupl)," rows are duplicated.\nAll rows in 'ldf' must be unique, when all columns except [",value.col,"] are considered.")
  ldf2 = ldf0 = ldf

  # If required, remove all rows with pre-existing NA in the value.col
  if(ignore.na) {
    ldf2 = ldf2[!is.na(ldf2[,value.col]),]
    if(notify & nrow(ldf) > nrow(ldf2)) message(nrow(ldf) - nrow(ldf2), " rows (",(nrow(ldf) - nrow(ldf2))/nrow(ldf),") were excluded before trimming, becaused they contained NA.")}

  # If required select or deselect specific rows
  ldf2 = apply.to.not.to(ldf = ldf2, apply.to = apply.to, apply.not.to = apply.not.to)

  # If warn.ties is enabled, issue a warning about tied values: Create a new df with the number and proportion of unique values in each grouping and if there are any ties (non-unique values) issue a warning
  if(warn.ties){
    ldf3 = ldf2[!is.na(ldf2[,value.col]),]
    unik = suppressMessages({ldf3 %>% dplyr::summarise(unique_values = length(unique(.data[[value.col]])), prop = length(unique(.data[[value.col]])) / dplyr::n())})
    non.unik = unik$unique_values < dplyr::group_size(ldf3)
    if(any(non.unik)) if(trim.ties) warning(sum(non.unik)," of ",length(non.unik)," groupings contain tied values. Due to trim.ties=T, all instances of a tied value crossing the trim-threshold are trimmed.") else warning(sum(non.unik)," of ",length(non.unik)," groupings contain tied values. Due to trim.ties=F, instances of tied values to be trimmed are sampled randomly. Can lead to unexpected behaviour when ties around the lower and upper trim-threshold are identical.")
    if(view.ties) return(unik[non.unik,])}

  # Shuffle the dataframe if edge ties should not be trimmed together; then slice min and max rows as defined by trim (either as n or prop) and store these rows in bot and top
  if(!trim.ties) ldf2 = ldf2 %>% dplyr::arrange(sample(dplyr::n())); if(trim[1]>=1 | trim[1]==0) { bot = ldf2 %>% dplyr::slice_min(n = trim[1], order_by = .data[[value.col]], with_ties = trim.ties); warn.b = F } else if (trim[1]>0) {bot = ldf2 %>% dplyr::slice_min(prop = trim[1], order_by = .data[[value.col]], with_ties = trim.ties); warn.b = ldf2 %>% dplyr::group_size() * trim[1] < 1} else stop("'trim' can not be negative")
  if(!trim.ties) ldf2 = ldf2 %>% dplyr::arrange(sample(dplyr::n())); if(trim[2]>=1 | trim[2]==0) { top = ldf2 %>% dplyr::slice_max(n = trim[2], order_by = .data[[value.col]], with_ties = trim.ties); warn.t = F } else if (trim[2]>0) {top = ldf2 %>% dplyr::slice_max(prop = trim[2], order_by = .data[[value.col]], with_ties = trim.ties); warn.t = ldf2 %>% dplyr::group_size() * trim[2] < 1} else stop("'trim' can not be negative")
  if(any(warn.b) | any(warn.t)) warning("The chosen trim-proportion for the lower|upper threshold (",trim[1],"|",trim[2],") equates to less than 1 value in ",sum(warn.b),"|",sum(warn.t)," of ",length(warn.b)," groupings.")

  # rbind bot and top, remove the value.col, remove duplicate rows (can happen with trim.ties = F), add a column trim containing 'trim', merge with ldf to make sure rows from bt are correctly assigned to rows from ldf, keep all rows of ldf (all.x=T) in the process -> thus rows in ldf without equivalent in top-bot get a NA in the trim-column, find these NA and use the !BOOL to index and NA the appropriate rows in ldf
  bt = rbind(bot, top)
  bt = bt[,!colnames(bt) %in% value.col]
  bt = cbind(bt[!duplicated(bt), ], trim = 'trim')
  trim.BOOL = !is.na(merge(ldf, bt, all.x = T)$trim)
  ldf[trim.BOOL, value.col] = NA

  # Notification of trimming
  report = trim.statistics(ldf = ldf, ldf0 = ldf0, report.by = report.by, apply.to = apply.to, apply.not.to = apply.not.to, print.stats = notify, value.col = value.col)

  # Return with or without stats
  if(return.report) return(list(ldf = as.data.frame(ldf), report = report)) else return(as.data.frame(ldf))
}


#' Trim values beyond a defined spread in groupings of an ldf
#'
#' Groups an ldf by certain columns, calculates a freely definable centrality and spread statistic for each grouping, and trims values that are more than x times smaller or larger than the spread by overwriting them with NA.
#' @param ldf A long format dataframe with a column (default: 'value') holding the values, a column (default: 'measure', e.g. reaction time) describing the value, a column (default: 'observation') linking different measures of the same observation together, and at least one further column being used as identifier-column (such as 'participant id'). Use example.trial.data() for a correctly formatted long dataframe. Tip: you can use dplyr::pivot_longer to transform a wide dataframe into a long dataframe.
#' @param value.col Name of the column holding the values in the ldf. Defaults to 'value'.
#' @param group.by Identifier-columns that are used to create groupings, to which the test is applied.
#' @param apply.to,apply.not.to A named list of vectors to (un)select rows. The names of the columns are specified by the element names, and the identifiers in this column to which the trimming should (not) be applied are specified by the vector. Default to NULL, which means no (un)selection is applied to any column.
#' @param center.func,spread.func Name of a function that can take a vector and output a scalar. These function are used to summarise the values of each grouping and find their 'central point' and spread, respectively. Defaults to 'mean' and 'sd', respectively. 'na.rm = T' is passed as additional argument.
#' @param trim Vector of one or two values. Multiplied with the spread and subtracted/added to the central tendency to get the lower and upper bound. A single value is duplicated. Defaults to c(2.5,2.5).
#' @param report.by Identifier-columns that are used to create groupings, for which statistics are printed when notify is TRUE. Defaults to NULL, which means one value for the entire ldf is returned.
#' @param return.report,return.stats Adds the trimming report/statistics of centrality and spread to the trimmed ldf in a list. Defaults to FALSE.
#' @param notify Provides progress information. Defaults to TRUE.
#' @importFrom dplyr '%>%'
#' @keywords trimming
#' @export
#' @examples
#' # Trimming each participants RTs that are more than 2 SD beyond the mean of this participant
#' ldf = trim.by_spread(ldf = example.trial.data(), group.by = 'id', apply.to = list(measure = 'RT'), trim = 2)
#' # Trimming each participants RTs that are beyond 2 SD, but only above, not below (using Inf in 'trim')
#' ldf = trim.by_spread(ldf = example.trial.data(), group.by = 'id', apply.to = list(measure = 'RT'), trim = c(Inf,2))
#' # Do the same, but now for each session separately and report trimming statistics it in more detail
#' ldf = trim.by_spread(ldf = example.trial.data(), group.by = c('id','session'), apply.to = list(measure = 'RT'), trim = c(Inf,2), report.by = c('session','condition'))
trim.values.by_spread = function(ldf, group.by, apply.to = NULL, apply.not.to = NULL, center.func = 'mean', spread.func = 'sd', trim = 2.5, report.by = NULL, return.report = F, return.stats = F, notify = T, value.col = 'value'){

  # Check if args match the ldf
  stp = check.arguments(); if(!is.null(stp)) stop(stp)

  # Ensure 'trim' has length 2
  if(length(trim) > 2) stop("'trim' can contain 1 or two 2 values")
  trim = rep(trim,2)[1:2]

  # Group the ldf and create a second ldf, in which the trimming is done
  ldf = ldf %>% dplyr::group_by(dplyr::across(dplyr::all_of(group.by)))
  dupl = ldf %>% dplyr::select(-!!dplyr::sym(value.col)) %>% dplyr::group_by_all() %>% dplyr::group_size()
  if(any(dupl > 1)) stop(sum(dupl-1)," of ",length(dupl)," rows in 'ldf' are duplicated.\n  All rows must be unique, when all columns except [",value.col,"] are considered.")
  ldf2 = ldf0 = ldf

  # If required select or deselect specific rows
  ldf2 = apply.to.not.to(ldf = ldf2, apply.to = apply.to, apply.not.to = apply.not.to)

  # In each grouping calculate the centrality and spread
  stats = suppressMessages({ldf2 %>% dplyr::summarise(center = do.call(match.fun(center.func), list(!!dplyr::sym(value.col), na.rm = T)),
                                                      spread = do.call(match.fun(spread.func), list(!!dplyr::sym(value.col), na.rm = T)))})
  # Join these stats back with the ldf
  ldf2 = suppressMessages({dplyr::left_join(ldf2, stats)})
  # Find values that are below or above the treshold
  ldf2$trim = ldf2[[value.col]] < ldf2$center - trim[1] * ldf2$spread | ldf2[[value.col]] > ldf2$center + trim[2] * ldf2$spread
  # Join back with the original ldf, in case any rows were removed by apply.to and apply.not.to. These rows become NA in column trim and need to be FALSE
  ldf = suppressMessages({dplyr::left_join(ldf, ldf2)}); ldf$trim[is.na(ldf$trim)] = F
  # Trim values based on the trim-column and remove now unnecessary columns center, spread, and trim
  ldf[ldf$trim,value.col] = NA
  ldf = ldf %>% dplyr::select(-center, -spread, - trim)

  # Notification of trimming
  report = trim.statistics(ldf = ldf, ldf0 = ldf0, report.by = report.by, apply.to = apply.to, apply.not.to = apply.not.to, print.stats = notify, value.col = value.col)

  # Return with or without report and stats
  rt = list(ldf = as.data.frame(ldf))
  if(return.report) rt = c(rt, report = report)
  if(return.stats) rt = c(rt, stats = stats)
  if(return.report | return.stats) return(rt) else return(as.data.frame(ldf))
}



#' Copy trimming from one identifier to another
#'
#' Copy the trimming-status of values (source-values) to other values (target-values), that differ only in a single identifier.
#' @param ldf A long format dataframe with one column holding the values and at least one column holding identifiers that can be used to group the values. Compare example.trial.data().
#' @param value.col Name of the column holding the values in the ldf. Defaults to 'value'.
#' @param group.by Identifier-columns that are used to create groupings, within which the trimming-status is copied from one row to another. Defaults to NULL, which means all columns besides 'value.col' and 'copy.by' are used.
#' @param copy.by Identifier-column used to specify rows within each grouping from which to which the trimming status should be copied.
#' @param copy.from,copy.to Identifiers in the copy.by-column that specify the source-rows ('copy.from') and target-rows ('copy.to') for copying. In each grouping, if any or all (defined by 'if.all') of the source-values are NA, all target values are trimmed.
#' @param apply.to,apply.not.to A named list of vectors to (un)select rows. The names of the columns are specified by the element names, and the identifiers in this column to which the trimming should (not) be applied are specified by the vector. Default to NULL, which means no (un)selection is applied to any column.
#' @param if.all Requires that all values defined by 'copy.from' are NA, in order to copy NAs to 'copy.to'. Defaults to FALSE, which means the target-values are trimmed when any source-value is NA.
#' @param report.by Identifier-columns that are used to create groupings, for which statistics are printed when notify is TRUE. Defaults to NULL, which means one value for the entire ldf is returned.
#' @param return.report Adds the trimming report to the trimmed ldf in a list. Defaults to FALSE.
#' @param notify Provides progress information. Defaults to TRUE.
#' @importFrom dplyr '%>%'
#' @keywords trimming
#' @export
#' @examples
#' # Trim the RT of example.trial.data
#' ldf = trim.by_interval(ldf = example.trial.data(), group.by = c('id', 'session', 'measure'), apply.to = list(measure = 'RT'), trim = 0.1)
#' # Copy the trimming within the measure-column. With copy.from and copy.to being NULL, all identifiers (here: RT and ACC) are tested, and when the value of any measure-identifier in the grouping is NA, all values in the grouping (note that a grouping here includes the trial!) are set to NA
#' ldf2 = copy.trimming(ldf = ldf, copy.by = 'measure')
#' # Copy the trimming across all participants, separate for session and trial. In other words: for ALL participants, trim the value of trial x in session x if any participant has an NA-value in trial x in session x.
#' ldf2 = copy.trimming(ldf = ldf, group.by = c('session','trial','measure'), copy.by = 'id', report.by = c('session','measure'))
#' # Without grouping by measure, ACC-values are affected as well
#' ldf2 = copy.trimming(ldf = ldf, group.by = c('session','trial'), copy.by = 'id', report.by = c('session','measure'))
#' # Do the same, but only if ALL participants have an NA-value for a given trial (Given the data in example.trial.data(), likely no copying was applied).
#' ldf2 = copy.trimming(ldf = ldf, group.by = c('session','trial','measure'), copy.by = 'id', report.by = c('session','measure'), if.all = T)
#' # Copy.by id, but do not specify anything else. Column 'condition' is automatically included in 'group.by', thus more groupings exist and the effect of a single NA is thus limited to fewer rows.
#' ldf2 = copy.trimming(ldf = ldf, group.by = c('session','trial','measure'), copy.by = 'id', report.by = c('session','measure','condition'), apply.not.to = list(measure = 'ACC'))
#' ldf2 = copy.trimming(ldf = ldf, copy.by = 'id', report.by = c('session','measure','condition'), apply.not.to = list(measure = 'ACC'))
copy.trimming = function(ldf, group.by = NULL, copy.by, copy.from = NULL, copy.to = NULL, apply.to = NULL, apply.not.to = NULL, if.all = F, notify = T, report.by = NULL, return.report = F, value.col = 'value'){

  # Check if args match the ldf
  stp = check.arguments(off.limit.colnames = 'rows.to.be.trimmed'); if(!is.null(stp)) stop(stp)

  # Get the copy identifiers
  copy = as.vector(unique(ldf[[copy.by]]))
  if(is.null(group.by)) {group.by = colnames(ldf)[!colnames(ldf) %in% c(copy.by, value.col)]; message("'group.by' is NULL. Using all columns except [",copy.by,"] and [",value.col,"] for grouping.")}
  if(notify & (is.null(copy.from) | is.null(copy.to))) message(if(is.null(copy.from)) "'copy.from'",if(is.null(copy.from) & is.null(copy.to)) ",",if(is.null(copy.to)) "'copy.to'"," is NULL",". Using all identifiers in [",copy.by,"]: ", paste(copy, collapse = ','))
  if(is.null(copy.from)) copy.from = copy
  if(is.null(copy.to  )) copy.to   = copy

  # Group the ldf and create a second ldf, in which the trimming is done; save ldf0 for trim.statistics() to compare
  ldf0 = ldf = ldf %>% dplyr::group_by(dplyr::across(dplyr::all_of(group.by)))

  # If required, select and deselect specific rows
  ldf2 = apply.to.not.to(ldf = ldf, apply.to = apply.to, apply.not.to = apply.not.to)

  # Rows may not be duplicated
  # dupl = duplicated(apply(ldf2[,c(group.by, copy.by)], 1, paste, collapse=''))
  # if(any(dupl)) {print(ldf2[dupl,c(group.by, copy.by)]); stop("These rows are duplicated. Values must be uniquely identified by identifiers in 'group.by' and 'copy.by'")}
  # At least one grouping should
  dupl = duplicated(apply(ldf2[,group.by], 1, paste, collapse=''))
  if(!any(dupl)) {warning("No groupings as defined by 'group.by' contain more than 1 value. Thus, no action was performed.")}

  # Filter rows defined by copy.from and determine whether all/any of the values in each grouping are NA
  copyna = ldf2[ldf2[[copy.by]] %in% copy.from,]
  if(if.all) copyna = suppressMessages({copyna %>% dplyr::summarise(rows.to.be.trimmed = all(is.na(.data[[value.col]])),.groups = 'keep')}) else copyna = suppressMessages({copyna %>% dplyr::summarise(rows.to.be.trimmed = any(is.na(.data[[value.col]])),.groups = 'keep')})
  # Merge this trimming template with the original ldf; all rows that don't exist in copyna are made FALSE; i.e. not to be trimmed
  copyna = suppressMessages({dplyr::left_join(ldf, copyna)})
  copyna$rows.to.be.trimmed[is.na(copyna$rows.to.be.trimmed)] = F
  # NA the value.col in the rows specified by copy.to in the groupings where NA were found
  ldf[copyna$rows.to.be.trimmed & ldf2[[copy.by]] %in% copy.to, value.col] = NA

  # Notification of trimming
  report = trim.statistics(ldf = ldf, ldf0 = ldf0, report.by = report.by, apply.to = apply.to, apply.not.to = apply.not.to, print.stats = notify, value.col = value.col)

  # Return with or without report
  if(return.report) return(list(ldf = as.data.frame(ldf), report = report)) else return(as.data.frame(ldf))
}


#' Return statistics about NA-values per grouping in an ldf
#'
#' Calculate statistics about NA-values in a specified groupings of a long format dataframe, as well as optionally the difference to a prior state of the ldf (ldf0).
#' @param ldf A long format dataframe with one column holding the values and at least one column holding identifiers that can be used to group the values. Compare example.trial.data().
#' @param ldf0 A previous state of ldf, to calculate the difference in NA-values
#' @param value.col Name of the column holding the values in ldf. Defaults to 'value'.
#' @param report.by Identifier-columns that are used to create groupings, for which statistics are returned.
#' @param apply.to,apply.not.to A named list of vectors to (un)select rows. The names of the columns are specified by the element names, and the identifiers in this column to which the trimming should (not) be applied are specified by the vector. Default to NULL, which means no (un)selection is applied to any column.
#' @param print.stats Prints the stats instead of returning them. Defaults to FALSE.
#' @importFrom dplyr '%>%'
#' @keywords keyword
#' @export
trim.statistics = function(ldf, ldf0 = NULL, report.by, apply.to = NULL, apply.not.to = NULL, print.stats = F, value.col = 'value'){

  ldf = apply.to.not.to(ldf = ldf, apply.to = apply.to, apply.not.to = apply.not.to)
  stat = suppressMessages({ldf %>% dplyr::group_by(dplyr::across(dplyr::all_of(report.by))) %>% dplyr::summarise(mean.NA = mean(is.na(.data[[value.col]])))})

  if(!is.null(ldf0)){
    ldf0 = apply.to.not.to(ldf = ldf0, apply.to = apply.to, apply.not.to = apply.not.to)
    if(any(dim(ldf0) != dim(ldf))) stop("'ldf' and 'ldf0' must share the same dimensions")
    if(any(ldf != ldf0, na.rm = T)) stop("Identifiers or values in 'ldf' and 'ldf0' differ. They may only contain different numbers of NA")
    stat0 = suppressMessages({ldf0 %>% dplyr::group_by(dplyr::across(dplyr::all_of(report.by))) %>% dplyr::summarise(mean.NA.before = mean(is.na(.data[[value.col]])))})
    stat$mean.NA.before = stat0$mean.NA.before
    stat$mean.NA.trimmed = stat$mean.NA - stat$mean.NA.before
    stat = cbind(stat[!colnames(stat)=='mean.NA'], mean.NA = stat$mean.NA)}

  if(print.stats){
    message("Proportion of trimmed values",if(!is.null(apply.to) | !is.null(apply.not.to)) ",",if(!is.null(apply.to)) " only including selected rows",if(!is.null(apply.to) & !is.null(apply.not.to)) " and",if(!is.null(apply.not.to)) " excluding deselected rows",":")
    report = as.data.frame(t(stat))
    colnames(report) = NULL
    print(report)}

  return(as.data.frame(stat))
}



#' Visualize the data
#'
#' Takes an ldf, groups it, and plots the values in observation and returns number and proportion of trimmed values for each grouping for each measure
#' @param ldf A long format dataframe with one column holding the values and at least one column holding identifiers that can be used to group the values. Compare example.trial.data().
#' @param value.col Name of the column holding the values in the ldf. Defaults to 'value'.
#' @param group.by Identifier-columns that are used to create groupings that are plotted.
#' @param view.together,view.separate Identifier-columns from 'group.by' that are shown in a single plot vs. those that are shown in separate plots. Only one of the two can be defined. Defaults to NULL, in which case 'view.together' defaults to 'group.by' and a single plot is generated.
#' @param observation Identifier-column that is NOT in 'group.by' and identifies the values within each grouping. Defaults to NULL, in which case an observation-column is automatically generated.
#' @param sort.observation Sorts observations within each grouping independently. Only relevant for point and line plots. Ticks on the x-axis will be removed. Defaults to FALSE.
#' @param by.observation Flips grouping and observation, such that observations become groupings and groupings become observations. Defaults to FALSE.
#' @param type Implemented plot types are 'point', 'line', 'point_line', 'hist', 'density', 'violin', and 'box'. Defaults to 'point_line'.
#' @param add.names Show the grouping-identifiers next to one datapoint of the grouping in point, line and density plots. Can be 'top', 'bottom', or 'random'. Defaults to NULL, which means names are not shown.
#' @param group.sep Separator used when collapsing the identifiers of multiple columns defined by 'view.together'. Affects the legend, the x-axis ticks for violinplot and boxplot, and 'add.names'. Defaults to '.'.
#' @param label Label for the x-axis for density plot and histogram and for the y-axis for all others. Defaults to ''.
#' @param fill Fill groupings with color. Not relevant for point and line plots. Defaults to FALSE.
#' @param bins Bin-number for histogram. Defaults to 30.
#' @param legend Adds a legend. Defaults to FALSE.
#' @importFrom dplyr '%>%'
#' @keywords statistics, plotting, visualization
#' @export
#' @examples
#' # View each participant's RT for each trial. Plot each session separately
#' view.data(ldf = example.trial.data(seed = 1), group.by = c('id','session','measure'), view = list(measure = 'RT'), view.separate = 'session')
#' # Values in each grouping were numbered automatically. To use the actual identifiers in the 'trial'-column, use observation (as trials are already numbered 1:10 and sorted accordingly, only the axis label changes in this case)
#' view.data(ldf = example.trial.data(seed = 1), group.by = c('id','session','measure'), view = list(measure = 'RT'), view.separate = 'session', observation = 'trial')
#' # Sort the observations by value. Axis ticks therefore disappear, because points on the same point on the x-axis no longer necessarily are from the same observation.
#' view.data(ldf = example.trial.data(seed = 1), group.by = c('id','session','measure'), view = list(measure = 'RT'), view.separate = 'session', observation = 'trial', sort.observation = T)
#' # View the RT for each condition in each session together, i.e. plot each participant separately
#' view.data(ldf = example.trial.data(), group.by = c('id','session','condition'), view = list(measure = 'RT'), view.together = c('session','condition'))
#' # Use a histogram plot to show the distribution of RTs for each combination of participant and condition
#' view.data(ldf = example.trial.data(), group.by = c('id','condition'), view = list(measure = 'RT'), type = 'hist', fill = T)
#' # Do the same with a density plot and add each groupings name to the highest point for easier identification
#' view.data(ldf = example.trial.data(), group.by = c('id','condition'), view = list(measure = 'RT'), type = 'density', add.names = 'top')
#' # Same, but as violin plot or boxplot. For these two, grouping names are automatically present in the x-axis.
#' view.data(ldf = example.trial.data(), group.by = c('id','condition'), view = list(measure = 'RT'), type = 'violin')
#' view.data(ldf = example.trial.data(), group.by = c('id','condition'), view = list(measure = 'RT'), type = 'box')
#' # Compare the plots when the grouping and observation are flipped (by.observation = T), such that one becomes the other. Note, that trials are pooled in a particular way: 1 is not a participant's first trial, but a participant's first A trial and B trial, respectively. May or may not be meaningful.
#' view.data(ldf = example.trial.data(), group.by = c('id','condition'), view = list(measure = 'RT'), type = 'violin')
#' view.data(ldf = example.trial.data(), group.by = c('id','condition'), view = list(measure = 'RT'), type = 'violin', by.observation = T)
view.data = function(ldf, group.by, view = NULL, view.not = NULL, view.together = NULL, view.separate = NULL, observation = NULL, by.observation = F, sort.observation = F, type = 'point_line', value.col = 'value', add.names = NULL, group.sep = '.', label = '', fill = F, bins = 30, legend = F){

  # Check whether the dataframe and the args have the right prorperties and fit together
  stp = check.arguments(off.limit.colnames = c('single_grouping.column','grouping_observation.column')); if(!is.null(stp)) stop(stp)

  if(!is.null(view.together) & !is.null(view.separate)) stop("'view.together' and 'view.separate' are complementary and only one can be defined")
  if(any(!view.together %in% group.by)) stop("'Only columns defined by 'group.by' can be used in 'view.together''")
  if(any(!view.separate %in% group.by)) stop("'Only columns defined by 'group.by' can be used in 'view.separate''")
  if(is.null(view.together) & is.null(view.separate)) view.together = group.by else if(is.null(view.together)) view.together = group.by[!group.by %in% view.separate] else if(length(view.together) < length(group.by)) view.separate = group.by[!group.by %in% view.together]
  if(!is.numeric(ldf[[value.col]])) warning("[",value.col,"] is not numeric. Can lead to unexpected behaviour for sorting and plotting.")

  # Extract relevant columns and group
  ldf2 = apply.to.not.to(ldf = ldf, apply.to = view, apply.not.to = view.not)
  # ldf2 = ldf2[,c(group.by, value.col)] %>% dplyr::group_by(dplyr::across(dplyr::all_of(group.by)))
  ldf2 = ldf2 %>% dplyr::group_by(dplyr::across(dplyr::all_of(group.by)))
  # If sorting is enabled for point/line plots, sort the ldf by value
  if(!by.observation & sort.observation & type %in% c('point','line','point_line')) ldf2 = ldf2 %>% dplyr::arrange(!!dplyr::sym(value.col))
  # Create a new column with a number for each row in the grouping
  ldf2 = ldf2 %>% dplyr::mutate(grouping_observation.column = dplyr::row_number())
  # If observation is defined, sorting is disabled, and type is point or line, replace identifiers in grouping_observation.column with the identifiers from observation
  if(!is.null(observation) & !sort.observation & type %in% c('point','line','point_line')) ldf2$grouping_observation.column = ldf2[[observation]]
  # Collapse identifiers in grouping cols to be viewed together into a single col
  ldf2$single_grouping.column = interaction(ldf2[,view.together], sep = group.sep)
  # If view.separate is defined, split the ldf. Store in a list to loop through
  if(!is.null(view.separate)) ldf2.list = split(ldf2, as.list(ldf2[,view.separate])) else ldf2.list = list('1' = ldf2)
  plots = list()
  # Loop through the list, make a plot, and store it in the plots-list
  for (i in 1:length(ldf2.list)) {

    # Create a ggplot and then fill it with data
    g = ggplot2::ggplot()

    if(!by.observation){
      # All plots except point and line plots need levels in the grouping (here grouping_observation.column) to plot properly
      if(type %in% c('box','violin','hist','density')) ldf2.list[[i]]$single_grouping.column = factor(ldf2.list[[i]]$single_grouping.column, levels = unique(ldf2.list[[i]]$single_grouping.column))

      # Show observation on x-axis, color by single_grouping.column
      if(any(type %in% c('point','point_line'))) g = g + ggplot2::geom_point(data = ldf2.list[[i]], ggplot2::aes(x = grouping_observation.column, y = !!dplyr::sym(value.col), color = single_grouping.column), show.legend = if(legend) T else F)
      if(any(type %in% c('line','point_line')))  g = g + ggplot2::geom_line( data = ldf2.list[[i]], ggplot2::aes(x = grouping_observation.column, y = !!dplyr::sym(value.col), color = single_grouping.column), show.legend = if(legend) T else F)

      # Show single_grouping.column on x-axis, color by single_grouping.column
      if(any(type == 'box'))    g = g + ggplot2::geom_boxplot(data = as.data.frame(ldf2.list[[i]]), ggplot2::aes(x = single_grouping.column, y = !!dplyr::sym(value.col), color = single_grouping.column, fill = if(fill) single_grouping.column), show.legend = if(legend) T else F)
      if(any(type == 'violin')) g = g + ggplot2::geom_violin( data = as.data.frame(ldf2.list[[i]]), ggplot2::aes(x = single_grouping.column, y = !!dplyr::sym(value.col), color = single_grouping.column, fill = if(fill) single_grouping.column), show.legend = if(legend) T else F)

      # Show measure on x-axis, color by single_grouping.column
      if(type == 'hist')      g = g + ggplot2::geom_histogram(data = ldf2.list[[i]], ggplot2::aes(x = !!dplyr::sym(value.col), color = single_grouping.column, fill = if(fill) single_grouping.column), bins = bins, show.legend = if(legend) T else F)
      if(type == 'density')   g = g + ggplot2::geom_density(  data = ldf2.list[[i]], ggplot2::aes(x = !!dplyr::sym(value.col), color = single_grouping.column, fill = if(fill) single_grouping.column), alpha = 0.1, show.legend = if(legend) T else F)
    }else{
      # All plots except point and line plots need levels in the grouping (here grouping_observation.column) to plot properly
      if(type %in% c('box','violin','hist','density')) ldf2.list[[i]]$grouping_observation.column = factor(ldf2.list[[i]]$grouping_observation.column, levels = unique(ldf2.list[[i]]$grouping_observation.column))

      # Show single_grouping.column on x-axis, color by observation
      if(any(type %in% c('point','point_line'))) g = g + ggplot2::geom_point(data = ldf2.list[[i]], ggplot2::aes(x = single_grouping.column, y = !!dplyr::sym(value.col), group = grouping_observation.column, color = as.factor(grouping_observation.column)), show.legend = if(legend) T else F)
      if(any(type %in% c('line','point_line'))) g = g + ggplot2::geom_line(data = ldf2.list[[i]], ggplot2::aes(x = single_grouping.column, y = !!dplyr::sym(value.col), group = grouping_observation.column, color = as.factor(grouping_observation.column)), show.legend = if(legend) T else F)

      # Show observation on x-axis, color by observation
      if(any(type == 'box'))    g = g + ggplot2::geom_boxplot(data = as.data.frame(ldf2.list[[i]]), ggplot2::aes(x = grouping_observation.column, y = !!dplyr::sym(value.col), color = grouping_observation.column, fill = if(fill) grouping_observation.column), show.legend = if(legend) T else F)
      if(any(type == 'violin')) g = g + ggplot2::geom_violin( data = as.data.frame(ldf2.list[[i]]), ggplot2::aes(x = grouping_observation.column, y = !!dplyr::sym(value.col), color = grouping_observation.column, fill = if(fill) grouping_observation.column), show.legend = if(legend) T else F)

      # Show measure on x-axis, color by observation
      if(type == 'hist')      g = g + ggplot2::geom_histogram(data = ldf2.list[[i]], ggplot2::aes(x = !!dplyr::sym(value.col), color = grouping_observation.column, fill = if(fill) grouping_observation.column), bins = bins, show.legend = if(legend) T else F)
      if(type == 'density')   g = g + ggplot2::geom_density(  data = ldf2.list[[i]], ggplot2::aes(x = !!dplyr::sym(value.col), color = grouping_observation.column, fill = if(fill) grouping_observation.column), alpha = 0.1, show.legend = if(legend) T else F)}

    # Put the name on one point of each grouping by extracting one point (max, min, or random) from the data of the plot itself (contains the x value); the resulting df is ordered by the levels of single_grouping.column, add a col with the actual names, and then add them as text on the spots defined by x and y values
    if(!is.null(add.names) & any(type %in% c('density', 'point', 'line', 'point_line'))){
      plotdata = ggplot2::ggplot_build(g)$data[[1]] %>% dplyr::group_by(group)
      if('top' %in% add.names)    plotdata = plotdata %>% dplyr::slice_max(order_by = y)
      else if('bottom' %in% add.names) plotdata = plotdata %>% dplyr::slice_min(order_by = y)
      else if('random' %in% add.names) plotdata = plotdata %>% dplyr::slice_sample()
      else (stop("'add.names' must be be 'top', 'bottom', or 'random'"))
      plotdata$name = levels(factor(ldf2.list[[i]]$single_grouping.column)) # Add a column with names
      g = g + ggplot2::geom_text(data = plotdata, ggplot2::aes(x = x, y = y, label = name), size = 3)}


    # Labs, legend, title; differs between by.observation
    if(!by.observation){
      if(type %in% c('point','line','point_line')) g = g + ggplot2::labs(x = if(!is.null(observation)) observation, y = label, title = paste(paste(view.separate, ldf2.list[[i]][1,c(view.separate)]), collapse = ' | '), color = paste(view.together, collapse = '\n'), fill = if(fill) paste(view.together, collapse = '\n'))
      if(type %in% c('violin','box'))              g = g + ggplot2::labs(x = paste(view.together, collapse = ' '),  y = label, title = paste(paste(view.separate, ldf2.list[[i]][1,c(view.separate)]), collapse = ' | '), color = paste(view.together, collapse = '\n'), fill = if(fill) paste(view.together, collapse = '\n'))
      if(type %in% c('density','hist'))            g = g + ggplot2::labs(x = label,                                            title = paste(paste(view.separate, ldf2.list[[i]][1,c(view.separate)]), collapse = ' | '), color = paste(view.together, collapse = '\n'), fill = if(fill) paste(view.together, collapse = '\n'))
    } else {
      if(type %in% c('point','line','point_line')) g = g + ggplot2::labs(x = paste(view.together, collapse = '\n'), y = label, title = paste(paste(view.separate, ldf2.list[[i]][1,c(view.separate)]), collapse = ' | '), color = observation, fill = if(fill) observation)
      if(type %in% c('violin','box'))              g = g + ggplot2::labs(x = observation,                           y = label, title = paste(paste(view.separate, ldf2.list[[i]][1,c(view.separate)]), collapse = ' | '), color = observation, fill = if(fill) observation)
      if(type %in% c('density','hist'))            g = g + ggplot2::labs(x = label,                                            title = paste(paste(view.separate, ldf2.list[[i]][1,c(view.separate)]), collapse = ' | '), color = observation, fill = if(fill) observation)}

    # When sorting was applied, the scale for point and line can not have ticks, because each position means something different
    if(sort.observation & type %in% c('point','line','point_line')) g = g + ggplot2::scale_x_discrete(labels = NULL, breaks = NULL)
    if(fill) g = g + ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(color = NULL)))
    if(!legend) g = g + ggplot2::theme(legend.position = 'none')
    g = g + ggplot2::theme_minimal()
    plots[[i]] = g}

  return(plots)}

