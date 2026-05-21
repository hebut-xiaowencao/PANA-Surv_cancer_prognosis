import numpy as np

# Concordance Index 
def concordance_index(y_true, y_pred):
    y_true = np.squeeze(y_true)
    y_pred = np.squeeze(y_pred)
    t = np.abs(y_true)
    e = (y_true > 0).astype(np.int32)

    concordant = 0.0
    permissible = 0.0
    n = len(t)

    for i in range(n):
        for j in range(i + 1, n):
            if t[i] == t[j]:
                continue
            if e[i] == 1 and t[i] < t[j]:
                permissible += 1
                if y_pred[i] < y_pred[j]:
                    concordant += 1
                elif y_pred[i] == y_pred[j]:
                    concordant += 0.5
            elif e[j] == 1 and t[j] < t[i]:
                permissible += 1
                if y_pred[j] < y_pred[i]:
                    concordant += 1
                elif y_pred[i] == y_pred[j]:
                    concordant += 0.5

    if permissible == 0:
        return np.nan
    return concordant / permissible
