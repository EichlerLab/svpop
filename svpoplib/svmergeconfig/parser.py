import ply.yacc

from . import lexer

class MergeParser(object):

    tokens = lexer.MergeLexer.tokens

    def __init__(self, **kwdargs):
        if 'write_tables' not in kwdargs.keys():
            kwdargs['write_tables'] = False

        if 'debug' not in kwdargs.keys():
            kwdargs['debug'] = False

        self.lexer = lexer.MergeLexer().lexer
        self.parser = ply.yacc.yacc(module=self, **kwdargs)

    def p_merge_config(self, p):
        """
        merge_config : KEYWORD
                     | KEYWORD ':' ':'
                     | KEYWORD ':' ':' spec_list
        """

        if len(p) == 2:
            spec_list = []
        elif len(p) == 5:
            spec_list = p[4]
        else:
            raise RuntimeError(f'BUG: p_merge_config: Expected 2 or 5 elements, found {len(p)}')

        p[0] = {
            'strategy': p[1],
            'spec_list': spec_list
        }

    # spec_list
    def p_spec_list(self, p):
        """
        spec_list : spec
                  | spec ':' spec_list
        """
        if len(p) == 2:
            p[0] = p[1]
        else:
            p[0] = p[1] + p[3]

    def p_spec(self, p):
        """
        spec : KEYWORD
             | KEYWORD '(' val_list ')'
             | t_match
        """

        is_match = \
            len(p) == 2 and \
            issubclass(p[1].__class__, list) and \
            len(p[1]) > 0 and \
            issubclass(p[1][0].__class__, tuple) and \
            len(p[1][0])> 0 and \
            p[1][0][0] == 'match'

        if is_match:
            p[0] = [{
                'type': 'match',
                'val_list': p[1][0][1]
            }]
        else:
            p[0] = [{
                'type': p[1],
                'val_list': p[3] if len(p) == 4 else []
            }]

    ### Value primitives ###
    # def p_val_primitive_list(self, p):
    #     """
    #     val_primitive_list : val_primitive
    #                        | val_primitive ',' val_primitive_list
    #     """
    #
    #     if len(p) == 2:
    #         p[0] = p[1]
    #     else:
    #         p[0] = p[1] + p[3]

    # spec_primitive and primitive types
    #
    # Tuples of:
    # 1) type
    # 2) value (None if unlimited)
    # 3) key (if key=val, else None)
    def p_val_primitive(self, p):
        """
        val_primitive : t_int
                      | t_float
                      | t_unlimited
                      | t_match
        """
        p[0] = p[1]

    def p_t_int(self, p):
        """
        t_int : T_INT
              | T_INT_MULT
        """
        p[0] = [('int', p[1], None)]

    def p_t_float(self, p):
        """
        t_float : T_FLOAT
                | T_FLOAT_EXP
        """
        p[0] = [('float', p[1], None)]

    def p_t_unlimited(self, p):
        """
        t_unlimited : T_UNLIMITED
        """
        p[0] = [('unlimited', None, None)]


    #
    # Values (keyed or primitive)
    #
    def p_val_list(self, p):
        """
        val_list : val
                 | val ',' val_list
        """

        if len(p) == 2:
            p[0] = p[1]
        else:
            p[0] = p[1] + p[3]

    # spec and spec_list
    def p_val(self, p):
        """
        val : val_primitive
            | KEYWORD '=' val_primitive
            | KW_MATCH '=' val_primitive
        """

        if len(p) == 2:
            p[0] = p[1]
        elif len(p) == 4:
            p[0] = [(p[3][0][0], p[3][0][1], p[1])]
        else:
            raise RuntimeError(f'Parser bug: Expected 2 or 4 elements for "spec" rule, found {len(p)}')

    # match
    def p_t_match(self, p):
        """
        t_match : KW_MATCH '(' val_list ')'
                | KW_MATCH
        """

        if len(p) == 5:
            match_list = p[3]
        elif len(p) == 2:
            match_list = []
        else:
            raise RuntimeError(f'Parser bug: Expected 2 or 5 elements for "match" rule, found {len(p)}')

        p[0] = [('match', match_list, None)]

    # Error handling
    def p_error(self, p):
        self.parser.errtok = p

        if p is not None:
            raise RuntimeError(
                'Syntax error at position {} ("{}"): {}'.format(
                    p.lexpos,
                    p.value,
                    '"' + (p.lexer.lexdata[p.lexpos:(p.lexpos + 20)] + '...' if len(p.lexer.lexdata) - p.lexpos > 20 else '') + '"'
                        if len(p.lexer.lexdata) - p.lexpos > 5 else 'at end of input'
                )
            )
        else:
            # p is None if an error occurs at the end
            raise RuntimeError(
                'Syntax error at end of format string (incomplete expression): Possibly missing a closing ")"?'
            )

    # # Build
    # def build(self, **kwdargs):
    #     if 'write_tables' not in kwdargs.keys():
    #         kwdargs['write_tables'] = False
    #
    #     if 'debug' not in kwdargs.keys():
    #         kwdargs['debug'] = False
    #
    #     self.lexer = lexer.MergeLexer().lexer
    #     self.parser = ply.yacc.yacc(module=self, **kwdargs)


    # Parse
    def parse(self, *args, **kwdargs):
        return self.parser.parse(*args, **kwdargs)
