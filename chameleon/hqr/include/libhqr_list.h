/**
 *
 * @file libhqr_list.h
 *
 * List module for the adapted reduction tree algorithm.
 *
 * @copyright 2017-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2018-09-30
 *
 */
#ifndef _libhqr_list_h_
#define _libhqr_list_h_

struct libhqr_list_elt_s;
typedef struct libhqr_list_elt_s  libhqr_list_elt_t;
typedef struct libhqr_list_elt_s *libhqr_list_t;

/**
 * @brief A list element data structure
 */
struct libhqr_list_elt_s {
    libhqr_list_elt_t *next; /**< The pointer to the next element in the list  */
    int id;                  /**< The id of the element                        */
    int date;                /**< The date of last modification of the element */
};

/**
 * @brief Push an element into a sorted list
 * @param[in,out] list_head
 *             The head of the list.
 * @param[in,out] elt
 *             The element to add to the list at the right place.
 */
static inline void
libhqr_list_push( libhqr_list_t *list_head, libhqr_list_elt_t *elt )
{
    libhqr_list_elt_t *prev, *next;
    assert( list_head != NULL );
    assert( (elt != NULL) && (elt->next == NULL) );

    prev = NULL;
    next = *list_head;
    while( (next != NULL) &&
           (next->date <= elt->date) )
    {
        prev = next;
        next = prev->next;
    }

    /* We replaced the head */
    if ( prev == NULL ) {
        *list_head = elt;
    }
    else {
        prev->next = elt;
    }
    elt->next = next;
}

/**
 * @brief Pop the first element of a sorted list
 * @param[in,out] list_head
 *             The head of the list. On exit, the head is removed.
 * @return The first element of the list.
 */
static inline libhqr_list_elt_t *
libhqr_list_pop( libhqr_list_t *list_head )
{
    libhqr_list_elt_t *elt;
    assert( list_head != NULL );
    elt = *list_head;
    if ( elt != NULL ) {
        *list_head = elt->next;
        elt->next = NULL;
    }
    return elt;
}

#endif /* _libhqr_list_h_ */
