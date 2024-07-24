# Program to print the velocity of each atom given two time-steps of a trajectory.xyz file
# 
# These trajectory.xyz files were generated via TDDFT simulations with the same Boltzmann distribution
# atom velocity method for their random velocities. 
# Program is useful for comparing initial velocities between what was originally generated
def read_positions(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    positions = []
    times = []
    start = 0
    while start < len(lines):
        # Read the number of atoms
        num_atoms = int(lines[start].strip())
        start += 1

        # Read the time step information
        time_info = lines[start].strip().split()
        time = float(time_info[-1])
        times.append(time)
        start += 1

        # Read the positions
        current_positions = []
        for _ in range(num_atoms):
            line = lines[start].strip().split()
            atom_type = line[0]
            x, y, z = float(line[1]), float(line[2]), float(line[3])
            current_positions.append((atom_type, x, y, z))
            start += 1

        positions.append(current_positions)

    return positions, times

def calculate_velocity(pos1, pos2, time1, time2):
    velocities = []
    for i in range(len(pos1)):
        _, x1, y1, z1 = pos1[i]
        _, x2, y2, z2 = pos2[i]
        vx = (x2 - x1) / (time2 - time1)
        vy = (y2 - y1) / (time2 - time1)
        vz = (z2 - z1) / (time2 - time1)
        velocities.append((vx, vy, vz))
    return velocities

def print_velocities(velocities):
    for _, (vx, vy, vz) in enumerate(velocities):
        print(f"{vx:.5e}\t{vy:.5e}\t{vz:.5e}")

def main():
    filename = 'positions.xyz'
    positions, times = read_positions(filename)

    if len(positions) != 2 or len(times) != 2:
        print("Error: The file must contain exactly two time steps.")
        return

    pos1, pos2 = positions
    time1, time2 = times

    velocities = calculate_velocity(pos1, pos2, time1, time2)
    print_velocities(velocities)

if __name__ == "__main__":
    main()
