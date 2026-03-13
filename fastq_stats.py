def main():
    filepath = 'reads.fastq'
    
    # Инициализация переменных для базовой статистики
    total_reads = 0
    total_length = 0
    min_length = float('inf')
    max_length = 0
    gc_count = 0
    
    # Статистика по качеству (позиция 10)
    total_quality_pos_10 = 0
    count_pos_10 = 0

    reads_removed_cnt = 0 # Количество удаленных прочтений (длина стала <= 0)
    filtered_lengths = [] # Длины оставшихся прочтений
    reads_filtered_len_60 = 0 # Оставшиеся с длиной >= 60
    
    # Параметры тримминга
    WINDOW_SIZE = 5
    QUAL_THRESHOLD = 30
    TOTAL_QUAL_LIMIT = WINDOW_SIZE * QUAL_THRESHOLD # Сумма 150
    
    try:
        with open(filepath, 'r') as f:
            line_idx = 0
            
            for line in f:
                line = line.strip()
                row_idx = line_idx % 4
                
                # Строка с последовательностью (индекс 1)
                if row_idx == 1:
                    seq_len = len(line)
                    
                    total_reads += 1
                    total_length += seq_len
                    gc_count += line.count('G') + line.count('C')
                    
                    if seq_len < min_length:
                        min_length = seq_len
                    if seq_len > max_length:
                        max_length = seq_len
                
                # Строка с качеством (индекс 3)
                elif row_idx == 3:
                    # Среднее качество на 10-й позиции (индекс 9)
                    if len(line) >= 10:
                        quality = ord(line[9]) - 33
                        total_quality_pos_10 += quality
                        count_pos_10 += 1
                    
                    quals = [ord(c) - 33 for c in line]
                    
                    kept_len = 0 # Если 0, значит прочтение отбрасывается
                    
                    # Проверяем, хватает ли длины для окна
                    if len(quals) >= WINDOW_SIZE:
                        # Сумма качества в первом окне
                        w_total = sum(quals[:WINDOW_SIZE])
                        
                        # Если первое окно хорошее (сумма > 150)
                        if w_total > TOTAL_QUAL_LIMIT:
                            kept_len = len(quals)
                            read_len = len(quals)
                            
                            # Проходим окном по прочтению
                            for i in range(read_len - WINDOW_SIZE):
                                # Обновляем сумму: вычитаем уходящий элемент, добавляем приходящий
                                w_total = w_total - quals[i] + quals[i+WINDOW_SIZE]
                                
                                # Если среднее качество <= 30 (сумма <= 150), обрезаем
                                if w_total <= TOTAL_QUAL_LIMIT:
                                    kept_len = i + WINDOW_SIZE
                                    break
                            
                            # Удаляем низкокачественные основания с конца (backtracking)
                            # Пока последний символ плохой (< 30) и длина > 1
                            while kept_len > 1 and quals[kept_len - 1] < QUAL_THRESHOLD:
                                kept_len -= 1
                            
                            # Если осталось меньше 1 основания, считаем прочтение удаленным
                            if kept_len < 1:
                                kept_len = 0
                        else:
                            # Первое окно плохое - удаляем прочтение
                            kept_len = 0
                    else:
                        # Слишком короткое прочтение - удаляем
                        kept_len = 0
                        
                    # Сбор статистики после тримминга
                    if kept_len == 0:
                        reads_removed_cnt += 1
                    else:
                        filtered_lengths.append(kept_len)
                        # Фильтрация по длине >= 60
                        if kept_len >= 60:
                            reads_filtered_len_60 += 1
    
                line_idx += 1
                
    except FileNotFoundError:
        print("Файл reads.fastq не найден")
        return

    # Вывод результатов
    if total_reads > 0:
        # Округление до ближайшего целого
        avg_len = int(total_length / total_reads + 0.5)
        gc_content = (gc_count / total_length) * 100
        
        avg_qual_10 = 0
        if count_pos_10 > 0:
            avg_qual_10 = int(total_quality_pos_10 / count_pos_10 + 0.5)
            
        # Статистика по отфильтрованным прочтениям
        min_filt = 0
        max_filt = 0
        avg_filt = 0
        if filtered_lengths:
            min_filt = min(filtered_lengths)
            max_filt = max(filtered_lengths)
            avg_filt = int(sum(filtered_lengths) / len(filtered_lengths) + 0.5)

        print(f"Общее число прочтений: {total_reads}")
        print(f"Минимальная длина прочтения: {min_length}")
        print(f"Средняя длина прочтения: {avg_len}")
        print(f"Максимальная длина прочтения: {max_length}")
        print(f"GC-состав: {gc_content:.2f}%")
        print(f"Среднее качество на позиции 10: {avg_qual_10}")
        
        print("-" * 20)
        print(f"Сколько прочтений подверглось триммингу (было удалено): {reads_removed_cnt}")
        print(f"Минимальная длина (в отфильтрованном файле): {min_filt}")
        print(f"Средняя длина (в отфильтрованном файле): {avg_filt}")
        print(f"Максимальная длина (в отфильтрованном файле): {max_filt}")
        print(f"Оставшееся число прочтений (с длиной >= 60): {reads_filtered_len_60}")

if __name__ == "__main__":
    main()
