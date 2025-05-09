import random

def generate_dna_sequence(length):
    return ''.join(random.choices('ACGT', k=length))

def insert_name(sequence, name):
    pos = random.randint(0, len(sequence))
    return sequence[:pos] + name + sequence[pos:]

def calculate_statistics(sequence):
    counts = {nuc: sequence.count(nuc) for nuc in 'ACGT'}
    total = sum(counts.values())
    percentages = {nuc: round((counts[nuc] / total) * 100, 1) for nuc in counts}
    cg_ratio = round(((counts['C'] + counts['G']) / (counts['A'] + counts['T'])) * 100, 1)
    return percentages, cg_ratio

def save_to_fasta(file_name, sequence_id, description, sequence):
    with open(file_name, 'w') as f:
        f.write(f">{sequence_id} {description}\n")
        f.write(sequence + "\n")

def main():
    length = int(input("Podaj długość sekwencji: "))
    sequence_id = input("Podaj ID sekwencji: ")
    description = input("Podaj opis sekwencji: ")
    name = input("Podaj imię: ")

    dna_sequence = generate_dna_sequence(length)
    sequence_with_name = insert_name(dna_sequence, name)

    save_to_fasta(f"{sequence_id}.fasta", sequence_id, description, sequence_with_name)

    percentages, cg_ratio = calculate_statistics(dna_sequence)

    print(f"Sekwencja została zapisana do pliku {sequence_id}.fasta")
    print("Statystyki sekwencji:")
    for nuc, perc in percentages.items():
        print(f"{nuc}: {perc}%")
    print(f"%CG: {cg_ratio}")

if __name__ == "__main__":
    main()
