import subprocess
import re
import csv


block_sizes = [8, 16, 32, 48, 50, 64, 72, 88, 100, 120]
iterations = [1, 10, 100, 1000, 5000, 10000, 50000, 100000]

pattern = r'The average time taken by (.*) method is (.*)ms for \d+ iterations\.'
outputs = []

for block_size in block_sizes:
    for iteration in iterations:
        matrix_size = block_size * 8
        command = [
            'bsub', '-I', '-b', '-q', 'q_ustc', '-N', '1', '-np', '1', '-cgsp', '64', '-cache_size', '0',
            './build/bin/stencil_main',
            '-m', str(matrix_size), '-b', str(block_size), '-i', str(iteration), '-w', '1'
        ]

        print(f'----------- block size: {block_size}, iteration: {iteration} -----------')
        result = subprocess.run(command, capture_output=True, text=True)
        output = result.stdout.strip()
        print(output)

        result_dict = {
            'Block Size': block_size,
            'Iteration': iteration,
        }

        matches = re.finditer(pattern, output)
        for match in matches:
            method, time = match.groups()
            result_dict[f'{method}'] = f"{float(time):.3f}"

        outputs.append(result_dict)


with open('output.csv', 'w', newline='') as file:
    fieldnames = outputs[0].keys()
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(outputs)
