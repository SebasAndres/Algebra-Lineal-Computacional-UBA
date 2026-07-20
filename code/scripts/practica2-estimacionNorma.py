import numpy as np

import matplotlib.pyplot as plt

N = 10
ITER = 100

def estimarNorma(A, iter=ITER):
	pred = 0
	time_series_pred = []
	for i in range(iter):
		x = np.random.random((N,1))
		pred = max(
			pred,
			np.linalg.norm(A@x) / np.linalg.norm(x)
		)
		time_series_pred.append(pred)
	return pred, time_series_pred

if __name__ == '__main__':

	A = np.random.random((N,N))
	real_norm = np.linalg.norm(A)
	pred_norm, time_series_pred = estimarNorma(A)

	print("Norma real: ", real_norm)
	print("Norma estimada: ", pred_norm)

	plt.plot(
		range(ITER),
		time_series_pred,
		color='blue'
	)

	plt.plot(
		range(ITER),
		[real_norm for _ in range(ITER)],
		color='red'
	)
	plt.show()
