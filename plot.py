import matplotlib.pyplot as plt

# Чтение данных из input.txt
x1 = []
y1 = []
with open('build/input.txt', 'r') as f:
    for line in f:
        data = line.split()
        x1.append(float(data[0]))
        y1.append(float(data[1]))

# Чтение данных из output.txt
x2 = []
y2 = []
with open('build/output.txt', 'r') as f:
    for line in f:
        data = line.split()
        x2.append(float(data[0]))
        y2.append(float(data[1]))

# Создание графика
plt.figure(figsize=(10, 6))

# График для исходных данных (input.txt)
plt.plot(x1, y1, label='Original Data (input.txt)', color='blue', marker='o')

# График для интерполированных данных (output.txt)
plt.plot(x2, y2, label='Interpolated Data (output.txt)', color='red', linestyle='--', marker='x')

# Настройка графика
plt.title('Interpolation of f(x) = ln(x) + e/x', fontsize=14)
plt.xlabel('x', fontsize=12)
plt.ylabel('f(x)', fontsize=12)
plt.legend()
plt.grid(True)

# Показать график
plt.savefig('plot.png')
