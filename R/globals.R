# variables exposed by RTMB::getAll() inside the objective (avoids R CMD check notes)
utils::globalVariables(c("logK","logr","logq","lbk","logm","logsdb","logsdi","tk","trho","logP","cdev","r.pr","K.pr","bk.pr",
  "q.pr","nq","q.uniform","est.m","m.pr","mu.m","est.sdb","sdb.ig","sdb.pr","spict","sdi.pr","dyn","est.k","kd.pr","catch.err",
  "TC","sdc","ebm","eb.rem","ebcc","stI","N","lagidx","biascor","rw","rw.scale","est.rho","rho.pr","rho.fix","P.bound","K.bound",
  "sdi.bound","nst","st.type","st.t","st.pr","eb.on","nI","obs_i","obs_t","ebci","sets.q","sets.var","W","SE2","xlim","logIobs"))
# ggplot2 aesthetics in the NUTS diagnostics
utils::globalVariables(c("iter","value","chain","bin","count","lag","acf","par","x","y","xmin","xmax","ymin","ymax","v","lt",
  "flag","what","density","metric","group","rk"))
