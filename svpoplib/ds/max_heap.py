"""
Max-heap implementation.

Heap operations occur in-place modifying the heap (list object). Elements inside the heap must be comparable with "<".
"""

import math


def _max_heapify(heap, i):
    """
    Create a max-heap from a specified index. Do not call directly, used by heapify().

    :param heap: Array to heapify.
    :param i: Heap from this index.
    """

    index_l = i * 2 + 1
    index_r = index_l + 1

    largest = i

    if index_l < len(heap) and heap[i] < heap[index_l]:
        largest = index_l

    if index_r < len(heap) and heap[largest] < heap[index_r]:
        largest = index_r

    if largest != i:
        heap[i], heap[largest] = heap[largest], heap[i]
        _max_heapify(heap, largest)

    return


def heapify(heap):
    """
    Create a max-heap from an unsorted list.

    :param heap: Array to heapify.
    """

    heap_len = len(heap)

    max_parent = math.floor((heap_len - 2) / 2)

    for i in range(max_parent, -1, -1):
        _max_heapify(heap, i)

    return


def pop(heap):
    """
    Remove and return the max item from a heap.

    :param heap: Heap to pop from.

    :return: Max value in heap

    :raises IndexError: If the heap is empty.
    """

    if len(heap) == 0:
        raise IndexError('Pop from empty heap')

    max_val = heap[0]

    heap[0] = heap[len(heap) - 1]
    heap.pop()

    i = 0

    while i < len(heap):
        index_l = i * 2 + 1
        index_r = index_l + 1

        largest = i

        if index_l < len(heap) and heap[i] < heap[index_l]:
            largest = index_l

        if index_r < len(heap) and heap[largest] < heap[index_r]:
            largest = index_r

        if largest != i:
            heap[i], heap[largest] = heap[largest], heap[i]
            i = largest
        else:
            i = len(heap)

    return max_val


def insert(heap, val):
    """
    Insert a value into the heap.

    :param heap: Heap list.
    :param val: Value to insert.
    """

    heap.append(val)

    i = len(heap) - 1
    i_parent = math.floor((i - 1) / 2)

    while i > 0 and heap[i_parent] < heap[i]:
        heap[i], heap[i_parent] = heap[i_parent], heap[i]
        i = i_parent
        i_parent = math.floor((i - 1) / 2)


def check_heap(heap):
    """
    Check heap structure. Throws an exception if any element contains a greater child element.

    :param heap: Heap.

    :raises RuntimeError: If any child elements are greater than their parent.
    """
    for i in range(len(heap)):
        index_l = i * 2 + 1
        index_r = index_l + 1

        if index_l >= len(heap):
            continue

        if index_l < len(heap) and heap[i] < heap[index_l]:
            raise RuntimeError(f'Value out of order at index parent {heap[i]} (index={i}) and child {heap[index_l]} (index={index_l})')

        if index_r < len(heap) and heap[i] < heap[index_r]:
            raise RuntimeError(f'Value out of order at index parent {heap[i]} (index={i}) and child {heap[index_r]} (index={index_r})')


def print_heap(heap):
    """
    Print a heap with the index before each item.

    :param heap: Heap to print.
    """
    for i in range(len(heap)):
        print(f'{i}: {heap[i]}')
