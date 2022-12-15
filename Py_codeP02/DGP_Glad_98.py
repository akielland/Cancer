# Simulating date with the DGP used in Glad 1998

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import AdaBoostRegressor


def true_model(x1, x2):
    # Generate y from
   return (2 + np.sin(2*np.pi*x1)) * (2+x2 - 2*x2**2 + 3*x2**5)

def noisy_model(x1, x2):
    # Generate y with added normal noise
    noise1 = np.random.normal(scale=0.5)
    noise2 = np.random.normal(scale=0.7)
    m1 = 2 + np.sin(2*np.pi*x1) + noise1
    m2 = 2 + x2 - 2*x2**2 + 3*x2**5 + noise2
    return m1 * m2

def DGP(size_x1, size_x2, noise=False):
    x1 = np.linspace(0, 1, size_x1)
    x2 = np.linspace(0, 1, size_x2)
    x1_x2_pairs = []
    y = []
    fun_vals = np.zeros([100, 100])
    for i, x_ in enumerate(x1):
        for j, y_ in enumerate(x2):
            if noise == False:
                fun_vals[i, j] = true_model(x_, y_)
                y.append(true_model(x_, y_))
            else:
                fun_vals[i, j] = noisy_model(x_, y_)
                y.append(noisy_model(x_, y_))
            x1_x2_pairs.append([x_, y_])
    X, Y = np.meshgrid(x1, x2)
    return x1_x2_pairs, y, fun_vals, X, Y


# linear model on x1 and x2
mechanistic = LinearRegression().fit(x1_x2_pairs, y_noisy)
y_predict = mechanistic.predict(x1_x2_pairs)
y_predict = y_predict.reshape(100, 100)

residuals = y_predict - fun_vals_noisy

clf = AdaBoostRegressor()
clf.fit(x1_x2_pairs, residuals.flatten())
adaboost = clf.fit(x1_x2_pairs, residuals.flatten())

ada_res = adaboost.predict(x1_x2_pairs)
full_prediction = y_predict - ada_res.reshape(100, 100)

# Linear model only on x1
xx1 = np.array(x1_x2_pairs)[:, 0]
mechanistic_2 = LinearRegression().fit(np.array([xx1]).transpose(), y_noisy)
y_predict_2 = mechanistic_2.predict(np.array([xx1]).transpose())
y_predict_2 = y_predict_2.reshape(100, 100)

residuals_2 = y_predict_2 - fun_vals_noisy

clf_2 = AdaBoostRegressor()
#clf_2.fit(x1_x2_pairs, residuals_2.flatten())
boost_2 = clf_2.fit(x1_x2_pairs, residuals.flatten())

res_boost_02 = boost_2.predict(x1_x2_pairs)
full_prediction_02 = y_predict - res_boost_02.reshape(100, 100)

## Boosting alone
clf_3 = AdaBoostRegressor()
boost_3 = clf_3.fit(x1_x2_pairs, y_noisy)
result_boost_03 = boost_3.predict(x1_x2_pairs)
result_boost_03 = result_boost_03.reshape(100, 100)


## Plot 1 ##
fig = plt.figure()
# set up the axes for the true model
ax = fig.add_subplot(2, 2, 1, projection='3d')
ax.set_zlim(0, 12)
x1_x2_pairs, y_true, fun_vals_true, X1, X2 = DGP(100, 100)
surf = ax.plot_surface(X1, X2, fun_vals_true, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# noisy simulation
ax = fig.add_subplot(2, 2, 2, projection='3d')
ax.set_zlim(0, 12)
x1_x2_pairs, y_noisy, fun_vals_noisy, X1, X2 = DGP(100, 100, noise=True)
surf = ax.plot_surface(X1, X2, fun_vals_noisy, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# prediction
ax = fig.add_subplot(2, 2, 3, projection='3d')
ax.set_zlim(0, 12)
surf = ax.plot_surface(X1, X2, full_prediction, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf = ax.plot_surface(X1, X2, y_predict, linewidth=0, color='black', alpha=0.2, antialiased=False)
# true model and prediction
ax = fig.add_subplot(2, 2, 4, projection='3d')
ax.set_zlim(0, 12)
surf = ax.plot_surface(X1, X2, fun_vals_true, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf = ax.plot_surface(X1, X2, full_prediction, linewidth=0, color='black', alpha=0.2, antialiased=False)


## Plot 2 ##
fig = plt.figure()
# set up the axes for the true model
ax = fig.add_subplot(2, 2, 1, projection='3d')
ax.set_zlim(0, 12)
x1_x2_pairs, y_true, fun_vals_true, X1, X2 = DGP(100, 100)
surf = ax.plot_surface(X1, X2, fun_vals_true, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# noisy simulation
ax = fig.add_subplot(2, 2, 2, projection='3d')
ax.set_zlim(0, 12)
x1_x2_pairs, y_noisy, fun_vals_noisy, X1, X2 = DGP(100, 100, noise=True)
surf = ax.plot_surface(X1, X2, fun_vals_noisy, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# prediction
ax = fig.add_subplot(2, 2, 3, projection='3d')
ax.set_zlim(0, 12)
surf = ax.plot_surface(X1, X2, full_prediction_02, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf = ax.plot_surface(X1, X2, y_predict_2, linewidth=0, color='black', alpha=0.2, antialiased=False)
# true model and prediction
ax = fig.add_subplot(2, 2, 4, projection='3d')
ax.set_zlim(0, 12)
surf = ax.plot_surface(X1, X2, fun_vals_true, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf = ax.plot_surface(X1, X2, full_prediction_02, linewidth=0, color='black', alpha=0.2, antialiased=False)

## plot 3 ##
fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.set_zlim(0, 12)
surf = ax.plot_surface(X1, X2, fun_vals_true, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf = ax.plot_surface(X1, X2, full_prediction_02, linewidth=0, color='black', alpha=0.2, antialiased=False)
ax = fig.add_subplot(1, 2, 2, projection='3d')
ax.set_zlim(0, 12)
surf = ax.plot_surface(X1, X2, fun_vals_true, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf = ax.plot_surface(X1, X2, result_boost_03, linewidth=0, color='black', alpha=0.2, antialiased=False)




surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
fig, ax = plt.subplots(221, subplot_kw={"projection": "3d"})



res = LinearRegression().fit(x1_x2_pairs, y)
y_predict = res.predict(x1_x2_pairs)
y_predict = y_predict.reshape(100, 100)

residuals = y_predict - noisy_simulation

clf = AdaBoostRegressor()
clf.fit(x1_x2_pairs, residuals.flatten())
adaboost = clf.fit(x1_x2_pairs, residuals.flatten())

adaboost.predict(x1_x2_pairs)

xx1 = np.array(x1_x2_pairs)[:, 0]
res2 = LinearRegression().fit(np.array([xx1]).transpose(), y)
y_predict2 = res2.predict(np.array([xx1]).transpose())