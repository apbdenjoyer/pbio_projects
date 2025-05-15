import os
# Importuje moduł os — pozwala na operacje systemowe, np. tworzenie folderów, przeglądanie katalogów itp.

import random
# i+Importuje moduł random — wykorzystywany do losowego generowania sekwencji i wybierania pozycji w niej

from pathlib import Path
# Importuje klasę Path z biblioteki pathlib — ułatwia pracę ze ścieżkami plików w sposób platformowo niezależny

import re
# Importuje moduł re — służy do pracy z wyrażeniami regularnymi (np. do wyszukiwania wzorców w sekwencjach DNA)


"""
READ ME:

[PL]

Rozszerzenie oryginalnego programu, stworzonego przez ChatGPT:
1. Dodano opcję zapisywania wielu sekwencji do jednego pliku
2. Dodano opcję wczytywania sekwencji z pliku (handling dla sekwencji o tym samym ID)
3. Dodano opcję (próby) odnalezienia zapisanego imienia/słowa w sekwencji 

Cel programu:
1. Możliwość generowania, zapisywania, odczytywania i generowania statystyk dla sekwencji DNA w formacie fasta

[EN]
Extension of the original program created by ChatGPT:
1. Added the option to save multiple sequences to a single file
2. Added the option to load sequences from a file (handling for sequences with the same ID)
3. Added the option to (attempt to) find a saved name/word within a sequence

Goal of program:
1. Ability to generate, save, write and get statistics of DNA sequences in fasta format
"""
# komentarz wieloliniowy zawierający opis rozszerzeń programu — zarówno po polsku jak i po angielsku


"""
ORIGINAL:
    brak/taki sam   /   none/same
MODIFIED: 
    [PL]
        (próba odnalezienia ukrytego słowa -> oczyszczenie sekwencji)
    [EN]
        (attempt at finding hidden word -> sequence cleaning)
"""
def find_name(sequence):
    # definiuje funkcję find_name, która próbuje znaleźć imię/słowo w sekwencji DNA

    """
        attempts to find a name.
        if the name contained either a/A, c/C, t/T, g/G then those will be counted as part of the original sequence.
    """

    match = re.search(r"^[ACTG]", sequence)
    # używa wyrażenia regularnego do znalezienia wzorca: od początku ciągu znajdź A, C, T lub G (jedną literę)

    return match[0], match.start(), re.sub(r"^[ACTG]", '', sequence)
    # Zwraca:
    # 1. znalezioną literę (nukleotyd) — match[0]
    # 2. pozycję tego znaku — match.start()
    # 3. sekwencję po usunięciu pierwszego znaku A/C/T/G — używa re.sub do zastąpienia tego znaku pustym ciągiem

"""
ORIGINAL:
    brak/taki sam   /   none/same
MODIFIED: 
"""
def generate_dna_sequence(length):
    # Generuje losową sekwencję DNA o zadanej długości
    # Używa funkcji random.choices do losowego wyboru nukleotydów (A, C, G, T)
    # Łączy wybrane nukleotydy w jeden ciąg znaków
    return ''.join(random.choices('ACGT', k=length))

"""
ORIGINAL:
    brak/taki sam   /   none/same
MODIFIED:
"""
def insert_name(sequence, name):
    # Wstawia podaną nazwę w losowe miejsce w sekwencji DNA
    pos = random.randint(0, len(sequence))  # losuje pozycję do wstawienia nazwy
    # Zwraca:
    # 1. zmodyfikowaną sekwencję (część przed pozycją + nazwa + część po pozycji)
    # 2. pozycję, w której wstawiono nazwę
    return sequence[:pos] + name + sequence[pos:], pos

"""
ORIGINAL:
    brak/taki sam   /   none/same
MODIFIED:
"""
def calculate_statistics(sequence):
    # Oblicza statystyki dla sekwencji DNA
    
    # Tworzy słownik z liczbą wystąpień każdego nukleotydu (A, C, G, T)
    counts = {nuc: sequence.count(nuc) for nuc in 'ACGT'}
    
    # Oblicza całkowitą długość sekwencji
    total = sum(counts.values())
    
    # Oblicza procentową zawartość każdego nukleotydu (zaokrąglone do 1 miejsca po przecinku)
    percentages = {nuc: round((counts[nuc] / total) * 100, 1) for nuc in counts}
    
    # Oblicza stosunek zawartości CG do AT (zaokrąglone do 1 miejsca po przecinku)
    cg_ratio = round(((counts['C'] + counts['G']) / (counts['A'] + counts['T'])) * 100, 1)
    
    # Zwraca:
    # 1. słownik z procentową zawartością każdego nukleotydu
    # 2. stosunek CG do AT
    return percentages, cg_ratio

"""
ORIGINAL:
    brak/taki sam   /   none/same
MODIFIED:
    [PL]
        (odczyt sekwencji z pliku -> by zapisane pliki mogły być aktualizowane o dodatkowe sekwencje)
    [EN]
        (sequence reading from file -> ability to append sequences to existing files)
"""
def read_fasta(file_name):
    # Wczytuje sekwencje z pliku w formacie FASTA
    with open(f"{FASTA_DIR}/{file_name}", 'r') as f:  # otwieramy plik do odczytu
        current_id = None  # aktualny identyfikator sekwencji
        current_description = None  # aktualny opis sekwencji
        current_sequence_lines = []  # lista linii aktualnej sekwencji

        for line in f:  # iterujemy po każdej linii pliku
            line = line.strip()  # usuwamy białe znaki z początku i końca
            if not line:  # pomijamy puste linie
                continue

            if line.startswith(">"):  # jeśli linia zaczyna się od ">" to jest to nagłówek sekwencji
                if current_id is not None:  # jeśli mamy poprzednią sekwencję, zapisujemy ją
                    full_sequence = ''.join(current_sequence_lines)  # łączymy wszystkie linie sekwencji
                    process_and_store_sequence(current_id, current_description, full_sequence)  # przetwarzamy i zapisujemy

                # przetwarzamy nowy nagłówek
                split = line[1:].split(' ', 1)  # dzielimy nagłówek na ID i opis
                current_id = split[0]  # pierwszy element to ID
                current_description = split[1] if len(split) > 1 else ""  # drugi (jeśli istnieje) to opis
                current_sequence_lines = []  # resetujemy listę linii sekwencji
            else:
                # jeśli to nie nagłówek, dodajemy linię do aktualnej sekwencji
                current_sequence_lines.append(line)

        # zapisujemy ostatnią sekwencję w pliku
        if current_id is not None:
            full_sequence = ''.join(current_sequence_lines)  # łączymy wszystkie linie
            process_and_store_sequence(current_id, current_description, full_sequence)  # przetwarzamy i zapisujemy

"""
ORIGINAL:
def save_to_fasta(file_name, sequence_id, description, sequence):
    os.makedirs("fasta_files", exist_ok=True)
    with open(f"fasta_files/{file_name}", 'w') as f:
        f.write(f">{sequence_id} {description}\n")
        f.write(sequence + "\n")
MODIFIED:
    [PL]
        (operowanie na kolekcji sekwencji -> możliwość zapisywania wielu sekwencji do jednego pliku)
    [EN]
        (using sequences collection -> ability to save multiple sequence per file)
"""
def write_fasta(file_name):
    # Zapisuje wszystkie sekwencje do pliku w formacie FASTA
    global sequences  # korzystamy z globalnej listy sekwencji
    
    with open(f"{FASTA_DIR}/{file_name}", 'w') as f:  # otwieramy plik do zapisu
        for seq in sequences:  # dla każdej sekwencji w liście
            # zapisujemy nagłówek (ID i opis)
            f.write(f">{seq['id']} {seq['description']}\n")
            # zapisujemy samą sekwencję
            f.write(seq['sequence_with_name']+"\n")


"""
ORIGINAL:
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
MODIFIED:
    [PL]
        (przeniesienie tworzenia sekwencji do create_sequence() -> modularyzacja)
        (przeniesienie wypisywania statystyk do print_statistics() -> modularyzacja)
        (dialog z użytkownikiem -> potrzeba przy rozbudowaniu aplikacji)
        (możliwość wczytywania sekwencji -> by zapisane pliki mogły być aktualizowane o dodatkowe sekwencje)
    [EN]
        (moved sequence creation to create_sequence() -> modularization)
        (moved statistics printing to print_statistics() -> modularization)
        (added user interaction -> necessary for expanding the application)
        (added the ability to load sequences -> so that saved files can be updated with additional sequences) 
"""
def print_statistics(seq):
    # Wyświetla statystyki dla pojedynczej sekwencji DNA
    percentages, cg_ratio = calculate_statistics(seq["sequence"])  # obliczamy procenty nukleotydów i stosunek CG
    output = [f"Sequence {seq['id']} statistics:",  # przygotowujemy nagłówek statystyk
              "Percentages:"]
    for nuc, perc in percentages.items():  # dla każdego nukleotydu
        output.append(f"\t*{nuc}: {perc}%")  # dodajemy jego procent
    output.append(f"\t*%CG: {cg_ratio}")  # dodajemy stosunek CG
    print("\n".join(output))  # wyświetlamy wszystkie statystyki

def create_sequence():
    # Tworzy nową sekwencję DNA na podstawie danych od użytkownika
    while True:
        try:
            length = int(input("Sequence length: "))  # pobieramy długość sekwencji
            break
        except ValueError:
            print("Err: Not a number.")  # obsługa błędu przy niepoprawnej wartości

    while True:
        sequence_id = input("Sequence ID: ")  # pobieramy identyfikator sekwencji
        
        # sprawdzamy czy ID jest unikalne
        unique = True
        for seq in sequences:
            if seq["id"] == sequence_id:
                unique=False
                break

        if not unique:
            print("Err: Sequence with this ID already exists.")
        else:
            break

    # pobieramy pozostałe dane
    description = input("Sequence description: ")  # opis sekwencji
    name = input("Name: ")  # nazwa do ukrycia w sekwencji
    dna_sequence = generate_dna_sequence(length)  # generujemy losową sekwencję DNA
    sequence_with_name, pos = insert_name(dna_sequence, name)  # wstawiamy nazwę do sekwencji
    return description, sequence_id, dna_sequence, sequence_with_name, name, pos

"""
ORIGINAL:
    brak/taki sam   /   none/same
MODIFIED:
    helper
"""
def choose_file_from_dir(directory):
    # Pomaga użytkownikowi wybrać plik z katalogu
    path = Path(directory)
    # tworzymy słownik z numerowanymi plikami
    found_files = {i: file for i, file in enumerate(os.listdir(path)) if os.path.isfile(path / file)}

    if not found_files:  # jeśli katalog jest pusty
        print(f"{directory} is empty. Aborting...")
        return None

    # wyświetlamy listę dostępnych plików
    print(f"Files inside directory '{directory}':")
    for i, file in found_files.items():
        print(f"\t*{file} ({i})")

    def is_valid_index(x):
        # sprawdza poprawność wybranego indeksu
        return x.isdigit() and int(x) in found_files

    # pobieramy wybór użytkownika
    index = int(get_valid_input("Choose a file (by number)\n$: ", is_valid_index, "Err: Invalid selection."))
    return found_files[index]

"""
ORIGINAL:
    brak/taki sam   /   none/same
MODIFIED:
    helper
"""
def get_valid_input(prompt, validator=lambda x: True, error_msg="Invalid input."):
    # Pobiera i sprawdza poprawność danych wejściowych
    while True:
        value = input(prompt)  # pobieramy dane
        if validator(value):  # sprawdzamy poprawność
            return value
        print(error_msg)  # wyświetlamy komunikat błędu

"""
ORIGINAL:
    brak/taki sam   /   none/same
MODIFIED:
    helper
"""
def confirm_choice(prompt):
    # Prosi użytkownika o potwierdzenie wyboru (Tak/Nie)
    while True:
        choice = input(f"{prompt} [Y/N]\n$: ").upper()  # pobieramy wybór
        if choice in ("Y", "N"):  # sprawdzamy poprawność
            return choice == "Y"
        print("Invalid input. Please enter Y or N.")


def main():
    # Główna funkcja programu
    global sequences
    os.makedirs(FASTA_DIR, exist_ok=True)  # tworzymy katalog na pliki FASTA

    while True:
        # wyświetlamy menu główne
        question = (f"Choose one of the following:\n"
                   f"\t*Exit program (1)"
                   f"\n\t*Load sequences from file (2)"
                   f"\n\t*Add new sequence to file (3)"
                   f"{'\n\t*Save file (4)' if len(sequences) > 0 else ''}\n$: ")

        choice = input(question)  # pobieramy wybór użytkownika
        match choice:
            case "1":  # wyjście z programu
                if sequences and not confirm_choice("Warning: unsaved sequences in memory. Are you sure you want to quit?"):
                    print("Cancelling...")
                    continue
                print("Quitting program...")
                exit(0)

            case "2":  # wczytywanie sekwencji z pliku
                file = choose_file_from_dir(FASTA_DIR)
                if file:
                    read_fasta(file)

            case "3":  # tworzenie nowej sekwencji
                description, sequence_id, sequence, sequence_with_name, name, pos = create_sequence()
                sequences.append({
                    "description": description,
                    "id": sequence_id,
                    "sequence": sequence,
                    "sequence_with_name": sequence_with_name,
                    "name": name,
                    "pos": pos
                })
                print(f"Info: Sequence {sequence_id} added.")

            case "4":  # zapisywanie do pliku
                filename = input("Please provide a filename (without .fasta)\n$: ")
                write_fasta(f"{filename}.fasta")
                print(f"Sequences saved to {filename}.fasta")
                
                for seq in sequences:  # wyświetlamy statystyki
                    print_statistics(seq)
                
                sequences = []  # czyścimy listę sekwencji

            case _:  # nieprawidłowy wybór
                print("Invalid command.")
            

"""
ORIGINAL:
    brak/taki sam   /   none/same
MODIFIED:
    helper
"""
def process_and_store_sequence(seq_id, description, sequence):
    # Przetwarza i zapisuje sekwencję w pamięci programu
    name, pos, sequence_clean = find_name(sequence)  # szukamy ukrytej nazwy

    # sprawdzamy czy wszystkie dane są poprawne
    if None in (seq_id, description, sequence, sequence_clean, name, pos):
        print("Error reading fasta")
        print((seq_id, description, sequence, sequence_clean, name, pos))
        return

    # tworzymy słownik z danymi sekwencji
    seq = {
        "description": description,
        "id": seq_id,
        "sequence": sequence_clean,
        "sequence_with_name": sequence,
        "name": name,
        "pos": pos
    }

    # sprawdzamy czy sekwencja o takim ID już istnieje
    if any(sq['id'] == seq_id for sq in sequences):
        while True:
            choice = input(f"Sequence {seq_id} already exists. Overwrite [Y] or Skip [N]: ").upper()
            if choice == "Y":  # nadpisujemy istniejącą sekwencję
                for i, sq in enumerate(sequences):
                    if sq['id'] == seq_id:
                        sequences[i] = seq
                        print(f"Sequence {seq_id} overwritten.")
                        print_statistics(seq)
                        break
                break
            elif choice == "N":  # pomijamy nową sekwencję
                print(f"Sequence {seq_id} skipped.")
                break
            else:
                print("Invalid choice. Please enter Y or N.")
    else:
        # dodajemy nową sekwencję
        sequences.append(seq)
        print(f"Sequence {seq_id} added.")
        print_statistics(seq)


# Zmienne globalne
sequences = []  # lista przechowująca wszystkie sekwencje
FASTA_DIR = Path(__file__).parent / "fasta_files"
# ścieżka katalogu z plikami FASTA

if __name__ == "__main__":
    main()  # uruchamiamy program