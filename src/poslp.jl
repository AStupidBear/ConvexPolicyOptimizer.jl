using JuQuant, OSQP, Elemental

scc_start()

nprocs() == 1 && using Mosek

pat = parseenv("pat", "bitmex_train")
prd = parseenv("prd", 1)

data = loaddata(pat * "*.h5")
data = downsample(data, prd)

savedata("data_train.h5", data)
w = fitadmm(nfeas(data), epochs = 10, α = 0.95, method = :LP)

p = Policy(model = ElRegressor(w = w), thresh = 0.5)
backtest(p, data, mode = "test")

scc_end()
