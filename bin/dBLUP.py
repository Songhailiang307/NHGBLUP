############################################################
# D2-GBLUP
#
# Strict DMU-style GBLUP
# + VanRaden2008 GRM
# + AI-REML
# + Henderson MME
# + Deep kernel extension
#
# METHODS:
#   tanh / mlp
#
# USAGE:
#   python d2blup.py tanh 0.2
#   python d2blup.py mlp 0.2
#
# OUTPUT:
#   Train_GEBV.txt
#   Validation_GEBV.txt
#   VanRaden2008_GRM.txt
#   Ghybrid.txt
#
# NOTE:
#   Output GRM matrices are generated for
#   ALL genotyped individuals
############################################################

import sys
import numpy as np
import pandas as pd

from scipy.linalg import solve
from scipy.optimize import minimize
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error
from sklearn.neural_network import MLPRegressor


############################################################
# 1. LOAD PLINK RAW
############################################################

def load_raw(raw_file):

    df = pd.read_csv(raw_file, sep=r"\s+")

    ids = df["IID"].astype(str).values

    geno = df.iloc[:, 6:].values.astype(np.float64)

    for j in range(geno.shape[1]):

        col = geno[:, j]

        mask = np.isnan(col)

        if np.sum(mask) > 0:

            mean_g = np.nanmean(col)

            col[mask] = mean_g

            geno[:, j] = col

    return geno, ids


############################################################
# 2. LOAD PHENOTYPE
############################################################

def load_pheno(pheno_file):

    pheno = pd.read_csv(pheno_file, sep=r"\s+", header=None)

    pheno.columns = ["ID", "PHENO"]

    pheno["ID"] = pheno["ID"].astype(str)

    return pheno


############################################################
# 3. BUILD DATASET
############################################################

def build_dataset(geno, geno_ids, pheno):

    geno_map = {iid: i for i, iid in enumerate(geno_ids)}

    pheno_ids = set(pheno["ID"])

    train_x, train_y, train_ids = [], [], []
    valid_x, valid_ids = [], []

    for _, row in pheno.iterrows():

        iid = row["ID"]

        if iid in geno_map:

            idx = geno_map[iid]

            train_x.append(geno[idx])

            train_y.append(float(row["PHENO"]))

            train_ids.append(iid)

    for iid in geno_ids:

        if iid not in pheno_ids:

            idx = geno_map[iid]

            valid_x.append(geno[idx])

            valid_ids.append(iid)

    return (
        np.array(train_x),
        np.array(train_y),
        train_ids,
        np.array(valid_x),
        valid_ids
    )


############################################################
# 4. VANRADEN GRM
############################################################

def vanraden_grm(X):

    p = np.mean(X, axis=0) / 2.0

    Z = X - 2.0 * p

    denom = np.sum(2.0 * p * (1.0 - p))

    if denom <= 0:
        denom = 1e-8

    G = (Z @ Z.T) / denom

    G += np.eye(G.shape[0]) * 1e-6

    return G, p, denom


############################################################
# 5. AI-REML
############################################################

class AIREML:

    def __init__(self, G, y):

        self.G = G
        self.y = y
        self.n = len(y)

    def reml_nll(self, params):

        vg = np.exp(params[0])
        ve = np.exp(params[1])

        V = vg * self.G + ve * np.eye(self.n)

        X = np.ones((self.n, 1))

        try:

            Vinv = np.linalg.inv(V)

            XtVinvX = X.T @ Vinv @ X

            beta = np.linalg.solve(
                XtVinvX,
                X.T @ Vinv @ self.y
            )

            r = self.y - X @ beta

            sign1, logdetV = np.linalg.slogdet(V)
            sign2, logdetX = np.linalg.slogdet(XtVinvX)

            nll = (
                0.5 * logdetV
                + 0.5 * logdetX
                + 0.5 * r.T @ Vinv @ r
            )

            return nll

        except:

            return 1e10

    def fit(self):

        init = np.log([1.0, 1.0])

        result = minimize(
            self.reml_nll,
            init,
            method="L-BFGS-B"
        )

        vg = np.exp(result.x[0])
        ve = np.exp(result.x[1])

        lam = ve / vg

        h2 = vg / (vg + ve)

        return vg, ve, h2, lam


############################################################
# 6. GBLUP
############################################################

class GBLUP_DMU:

    def __init__(self, G, lam):

        self.G = G
        self.Ginv = np.linalg.inv(G)
        self.lam = lam

    def fit(self, y):

        n = len(y)

        X = np.ones((n, 1))

        Z = np.eye(n)

        C11 = X.T @ X
        C12 = X.T @ Z
        C21 = Z.T @ X

        C22 = Z.T @ Z + self.lam * self.Ginv

        LHS = np.vstack([
            np.hstack([C11, C12]),
            np.hstack([C21, C22])
        ])

        RHS = np.concatenate([
            X.T @ y,
            Z.T @ y
        ])

        sol = solve(LHS, RHS)

        self.mu = sol[0]

        self.u = sol[1:]

        return self.u

    def predict_valid(self, G_valid):

        return G_valid @ self.Ginv @ self.u


############################################################
# 7. FEATURE TRANSFORM
############################################################

def feature_transform(train_x, valid_x, method, train_y=None):

    if method == "tanh":

        return np.tanh(train_x), np.tanh(valid_x)

    elif method == "mlp":

        mlp = MLPRegressor(
            hidden_layer_sizes=(128,),
            activation="tanh",
            solver="adam",
            max_iter=500,
            random_state=1
        )

        mlp.fit(train_x, train_y)

        W = mlp.coefs_[0]

        b = mlp.intercepts_[0]

        return (
            np.tanh(train_x @ W + b),
            np.tanh(valid_x @ W + b)
        )

    else:

        raise ValueError("method must be tanh or mlp")


############################################################
# 8. GRM FROM FEATURES
############################################################

def build_grm(E):

    G = (E @ E.T) / E.shape[1]

    G += np.eye(G.shape[0]) * 1e-6

    return G


############################################################
# 9. D2-GBLUP
############################################################

class D2GBLUP:

    def __init__(self, G_vr, G_deep, alpha, lam):

        self.G = (
            (1 - alpha) * G_vr
            + alpha * G_deep
        )

        self.Ginv = np.linalg.inv(self.G)

        self.alpha = alpha

        self.lam = lam

    def fit(self, y):

        n = len(y)

        X = np.ones((n, 1))

        Z = np.eye(n)

        C11 = X.T @ X
        C12 = X.T @ Z
        C21 = Z.T @ X

        C22 = Z.T @ Z + self.lam * self.Ginv

        LHS = np.vstack([
            np.hstack([C11, C12]),
            np.hstack([C21, C22])
        ])

        RHS = np.concatenate([
            X.T @ y,
            Z.T @ y
        ])

        sol = solve(LHS, RHS)

        self.mu = sol[0]

        self.u = sol[1:]

        return self.u

    def predict_valid(self, G_valid_vr, G_valid_deep):

        Gv = (
            (1 - self.alpha) * G_valid_vr
            + self.alpha * G_valid_deep
        )

        return Gv @ self.Ginv @ self.u


############################################################
# 10. MAIN
############################################################

if __name__ == "__main__":

    method = sys.argv[1]

    alpha = float(sys.argv[2])

    raw_file = "simu.raw"

    pheno_file = "pheno.txt"

    ########################################################
    # LOAD DATA
    ########################################################

    geno, geno_ids = load_raw(raw_file)

    pheno = load_pheno(pheno_file)

    train_x, train_y, train_ids, valid_x, valid_ids = build_dataset(
        geno,
        geno_ids,
        pheno
    )

    ########################################################
    # TRAINING GRM
    ########################################################

    G_vr, p, denom = vanraden_grm(train_x)

    ########################################################
    # FULL GRM (ALL INDIVIDUALS)
    ########################################################

    G_all_vr, _, _ = vanraden_grm(geno)

    ########################################################
    # AI-REML
    ########################################################

    aireml = AIREML(G_vr, train_y)

    vg, ve, h2, lam = aireml.fit()

    ########################################################
    # GBLUP
    ########################################################

    gblup = GBLUP_DMU(G_vr, lam)

    gblup_pred = gblup.fit(train_y)

    ########################################################
    # VALIDATION RELATIONSHIP
    ########################################################

    Z_train = train_x - 2 * p

    Z_valid = valid_x - 2 * p

    G_valid_vr = (
        Z_valid @ Z_train.T
    ) / denom

    valid_gblup = gblup.predict_valid(G_valid_vr)

    ########################################################
    # FEATURE TRANSFORM
    ########################################################

    E_train, E_valid = feature_transform(
        train_x,
        valid_x,
        method,
        train_y
    )

    ########################################################
    # DEEP GRM (TRAIN)
    ########################################################

    G_deep = build_grm(E_train)

    ########################################################
    # FULL DEEP GRM (ALL INDIVIDUALS)
    ########################################################

    if method == "tanh":

        E_all = np.tanh(geno)

    elif method == "mlp":

        mlp = MLPRegressor(
            hidden_layer_sizes=(128,),
            activation="tanh",
            solver="adam",
            max_iter=500,
            random_state=1
        )

        mlp.fit(train_x, train_y)

        W = mlp.coefs_[0]

        b = mlp.intercepts_[0]

        E_all = np.tanh(geno @ W + b)

    else:

        raise ValueError("method must be tanh or mlp")

    G_all_deep = build_grm(E_all)

    ########################################################
    # D2-GBLUP
    ########################################################

    d2 = D2GBLUP(
        G_vr,
        G_deep,
        alpha,
        lam
    )

    d2_pred = d2.fit(train_y)

    ########################################################
    # VALIDATION DEEP RELATIONSHIP
    ########################################################

    G_valid_deep = (
        E_valid @ E_train.T
    ) / E_train.shape[1]

    valid_d2 = d2.predict_valid(
        G_valid_vr,
        G_valid_deep
    )

    ########################################################
    # FULL HYBRID MATRIX
    ########################################################

    Ghybrid_all = (
        (1 - alpha) * G_all_vr
        + alpha * G_all_deep
    )

    ########################################################
    # SAVE MATRICES
    ########################################################

    np.savetxt(
        "VanRaden2008_GRM.txt",
        G_all_vr,
        fmt="%.6f"
    )

    np.savetxt(
        "Ghybrid.txt",
        Ghybrid_all,
        fmt="%.6f"
    )

    ########################################################
    # PRINT RESULTS
    ########################################################

    print(
        "\nGBLUP TRAIN:",
        pearsonr(train_y, gblup_pred)[0],
        mean_squared_error(train_y, gblup_pred)
    )

    print(
        "\nD2-GBLUP TRAIN:",
        pearsonr(train_y, d2_pred)[0],
        mean_squared_error(train_y, d2_pred)
    )

    ########################################################
    # SAVE TRAIN GEBV
    ########################################################

    pd.DataFrame({

        "ID": train_ids,
        "PHENO": train_y,
        "GBLUP": gblup_pred,
        "D2GBLUP": d2_pred

    }).to_csv(
        "Train_GEBV.txt",
        sep="\t",
        index=False
    )

    ########################################################
    # SAVE VALIDATION GEBV
    ########################################################

    pd.DataFrame({

        "ID": valid_ids,
        "GBLUP_GEBV": valid_gblup,
        "D2GBLUP_GEBV": valid_d2

    }).to_csv(
        "Validation_GEBV.txt",
        sep="\t",
        index=False
    )

    print("\nDone")