#include <config.h>

#include <dune/grid/io/file/dgfparser/dgfprojectionblock.hh>

namespace Dune
{

  namespace dgf
  {

    // ProjectionBlock
    // ---------------

    const char *ProjectionBlock::ID = "Projection";

    ProjectionBlock::ProjectionBlock ( std::istream &in, int dimworld )
    : BasicBlock( in, ID ),
      defaultFunction_( 0 )
    {
      while( getnextline() )
      {
        //std::cout << "Projection line:" << std::flush;
        nextToken();

        if( token.type == Token::functionKeyword )
        {
          nextToken();
          parseFunction();
        }
        else if( token.type == Token::defaultKeyword )
        {
          nextToken();
          parseDefault();
        }
        else if( token.type == Token::segmentKeyword )
        {
          nextToken();
          parseSegment();
        }
        else if( token.type != Token::endOfLine )
          DUNE_THROW( DGFException, "Error in " << *this << ": Invalid token." );
        matchToken( Token::endOfLine, "trailing tokens on line." );
      }
    }


    void ProjectionBlock::parseFunction ()
    {
      if( token.type != Token::string )
        DUNE_THROW( DGFException, "Error in " << *this << ": function name expected." );
      const std::string functionName = token.literal;
      if( functions_.find( functionName ) != functions_.end() )
        DUNE_THROW( DGFException, "Error in " << *this << ": redeclaration of function " << functionName << "." );
      nextToken();

      matchToken( Token::openingParen, "'(' expected." );
      if( token.type != Token::string )
        DUNE_THROW( DGFException, "Error in " << *this << ": variable name expected." );
      const std::string variableName = token.literal;
      nextToken();
      matchToken( Token::closingParen, "')' expected." );

      matchToken( Token::equals, "'=' expected." );
      const Expression *expression = parseExpression( variableName );

      //std::cout << std::endl << "Declaring function: " << functionName << "( " << variableName << " )" << std::endl;
      functions_[ functionName ] = expression;
    }


    const ProjectionBlock::Expression *
    ProjectionBlock::parseBasicExpression ( const std::string &variableName )
    {
      const Expression *expression;

      // parenthesized expression
      if( token.type == Token::openingParen )
      {
        nextToken();
        expression = parseExpression( variableName );
        matchToken( Token::closingParen, "')' expected." );
      }
      // vector constant
      else if( token.type == Token::openingBracket )
      {
        std::vector< double > value;
        nextToken();
        while( token.type == Token::number )
        {
          value.push_back( token.value );
          nextToken();
        }
        expression = new ConstantExpression( value );
        matchToken( Token::closingBracket, "']' expected." );
      }
      // norm expression
      else if( token.type == Token::normDelim )
      {
        nextToken();
        expression = new NormExpression( parseExpression( variableName ) );
        matchToken( Token::normDelim, "'|' expected." );
      }
      // number
      else if( token.type == Token::number )
      {
        expression = new ConstantExpression( token.value );
        nextToken();
      }
      // pi
      else if( token.type == Token::piKeyword )
      {
        expression = new ConstantExpression( M_PI );
        nextToken();
      }
      else if( token.type == Token::string )
      {
        if( token.literal != variableName )
        {
          FunctionMap::iterator it = functions_.find( token.literal );
          if( it == functions_.end() )
            DUNE_THROW( DGFException, "Error in " << *this << ": function " << token.literal << " not declared." );
          nextToken();
          matchToken( Token::openingParen, "'(' expected." );
          expression = new FunctionCallExpression( it->second, parseExpression( variableName ) );
          matchToken( Token::closingParen, "')' expected." );
        }
        else
        {
          expression = new VariableExpression;
          nextToken();
        }
      }
      else
        DUNE_THROW( DGFException, "Error in " << *this << ": " << "basic expression expected." );

      return expression;
    }


    const ProjectionBlock::Expression *
    ProjectionBlock::parsePostfixExpression ( const std::string &variableName )
    {
      const Expression *expression = parseBasicExpression( variableName );
      if( token.type == Token::openingBracket )
      {
        nextToken();
        if( (token.type != Token::number) || (double( int( token.value ) ) != token.value) )
          DUNE_THROW( DGFException, "Error in " << *this << ": integral number expected." );
        expression = new BracketExpression( expression, int( token.value ) );
        nextToken();
        matchToken( Token::closingBracket, "']' expected." );
      }
      return expression;
    }


    const ProjectionBlock::Expression *
    ProjectionBlock::parseUnaryExpression ( const std::string &variableName )
    {
      const Expression *expression;

      // unary minus
      if( (token.type == Token::additiveOperator) && (token.symbol == '-') )
      {
        nextToken();
        expression = new MinusExpression( parsePostfixExpression( variableName ) );
      }
      // sqrt
      else if( token.type == Token::sqrtKeyword )
      {
        nextToken();
        expression = new SqrtExpression( parseUnaryExpression( variableName ) );
      }
      // sin
      else if( token.type == Token::sinKeyword )
      {
        nextToken();
        expression = new SinExpression( parseUnaryExpression( variableName ) );
      }
      // cos
      else if( token.type == Token::cosKeyword )
      {
        nextToken();
        expression = new CosExpression( parseUnaryExpression( variableName ) );
      }
      else
        expression = parsePostfixExpression( variableName );

      return expression;
    }


    const ProjectionBlock::Expression *
    ProjectionBlock::parsePowerExpression ( const std::string &variableName )
    {
      const Expression *expression = parseUnaryExpression( variableName );
      while( token.type == Token::powerOperator )
      {
        nextToken();
        expression = new PowerExpression( expression, parseUnaryExpression( variableName ) );
      }
      return expression;
    }


    const ProjectionBlock::Expression *
    ProjectionBlock::parseMultiplicativeExpression ( const std::string &variableName )
    {
      const Expression *expression = parsePowerExpression( variableName );
      while( token.type == Token::multiplicativeOperator )
      {
        const char symbol = token.symbol;
        nextToken();
        if( symbol == '*' )
          expression = new ProductExpression( expression, parsePowerExpression( variableName ) );
        else if( symbol == '/' )
          expression = new QuotientExpression( expression, parsePowerExpression( variableName ) );
        else
          DUNE_THROW( DGFException, "Error in " << *this << ": Internal tokenizer error." );
      }
      return expression;
    }


    const ProjectionBlock::Expression *
    ProjectionBlock::parseExpression ( const std::string &variableName )
    {
      const Expression *expression = parseMultiplicativeExpression( variableName );
      while( token.type == Token::additiveOperator )
      {
        const char symbol = token.symbol;
        nextToken();
        if( symbol == '+' )
          expression = new SumExpression( expression, parseMultiplicativeExpression( variableName ) );
        else if( symbol == '-' )
          expression = new DifferenceExpression( expression, parseMultiplicativeExpression( variableName ) );
        else
          DUNE_THROW( DGFException, "Error in " << *this << ": Internal tokenizer error." );
      }
      return expression;
    }


    void ProjectionBlock::parseDefault ()
    {
      if( token.type != Token::string )
        DUNE_THROW( DGFException, "Error in " << *this << ": function name expected." );
      const std::string functionName = token.literal;
      nextToken();

      //std::cout << std::endl << "Default function: " << functionName << std::endl;
      FunctionMap::iterator it = functions_.find( functionName );
      if( it == functions_.end() )
        DUNE_THROW( DGFException, "Error in " << *this << ": function " << functionName << " not declared." );
      defaultFunction_ = it->second;
    }


    void ProjectionBlock::parseSegment ()
    {
      std::vector< unsigned int > faceId;
      while( token.type == Token::number )
      {
        if( double( (unsigned int)token.value ) != token.value )
          DUNE_THROW( DGFException, "Error in " << *this << ": integral number expected." );
        faceId.push_back( (unsigned int)token.value );
        nextToken();
      }
     
      if( token.type != Token::string )
        DUNE_THROW( DGFException, "Error in " << *this << ": function name expected." );
      const std::string functionName = token.literal;
      nextToken();

      //std::cout << std::endl << "Boundary projection for face";
      //for( size_t int i = 0; i < faceId.size(); ++i )
      //  std::cout << " " << faceId[ i ];
      //std::cout << ": " << functionName << std::endl;
      FunctionMap::iterator it = functions_.find( functionName );
      if( it == functions_.end() )
        DUNE_THROW( DGFException, "Error in " << *this << ": function " << functionName << " not declared." );
      boundaryFunctions_.push_back( std::make_pair( faceId, it->second ) );
    }


    void ProjectionBlock::matchToken ( const Token::Type &type, const std::string &message )
    {
      if( token.type != type )
        DUNE_THROW( DGFException, "Error in " << *this << ": " << message );
      if( type != Token::endOfLine )
        nextToken();
    }


    void ProjectionBlock::nextToken ()
    {
      int c;

      // eat white space
      while( ((c = line.peek()) == ' ') || (c == '\t') )
        line.get();

      // parse string literals
      if( ((c >= 'a') && (c <= 'z')) || ((c >= 'A') && (c <= 'Z')) )
      {
        token.type = Token::string;
        token.literal = "";
        while( ((c >= 'a') && (c <= 'z')) || ((c >= 'A') && (c <= 'Z')) )
        {
          token.literal += lowerCase( line.get() );
          c = line.peek();
        }

        if( token.literal == "default" )
          token.type = Token::defaultKeyword;
        else if( token.literal == "function" )
          token.type = Token::functionKeyword;
        else if( token.literal == "segment" )
          token.type = Token::segmentKeyword;
        else if( token.literal == "sqrt" )
          token.type = Token::sqrtKeyword;
        else if( token.literal == "sin" )
          token.type = Token::sinKeyword;
        else if( token.literal == "cos" )
          token.type = Token::cosKeyword;
        else if( token.literal == "pi" )
          token.type = Token::piKeyword;
      }
      // parse numeric constant
      else if( (c >= '0') && (c <= '9') )
      {
        token.type = Token::number;
        token.value = 0;
        while( (c >= '0') && (c <= '9') )
        {
          token.value = 10*token.value + double( c - '0' );
          token.literal += char( line.get() );
          c = line.peek();
        }
        if( c == '.' )
        {
          token.literal += line.get();
          c = line.peek();
          double factor = 0.1;
          while( (c >= '0') && (c <= '9') )
          {
            token.value += factor * double( c - '0' );
            token.literal += line.get();
            factor *= 0.1;
            c = line.peek();
          }
        }
      }
      // parse single character tokens
      else if( c == '=' )
        token.setSymbol( Token::equals, line.get() );
      else if( c == '(' )
        token.setSymbol( Token::openingParen, line.get() );
      else if( c == ')' )
        token.setSymbol( Token::closingParen, line.get() );
      else if( c == '[' )
        token.setSymbol( Token::openingBracket, line.get() );
      else if( c == ']' )
        token.setSymbol( Token::closingBracket, line.get() );
      else if( c == '|' )
        token.setSymbol( Token::normDelim, line.get() );
      else if( (c == '+') || (c == '-') )
        token.setSymbol( Token::additiveOperator, line.get() );
      else if( (c == '*') )
      {
        c = line.get();
        if( (line.peek() == '*') )
        {
          token.type = Token::powerOperator;
          line.get();
        }
        else
          token.setSymbol( Token::multiplicativeOperator, c );
      }
      else if( c == '/' )
        token.setSymbol( Token::multiplicativeOperator, line.get() );
      // parse end of line
      else if( c == std::stringstream::traits_type::eof() )
        token.type = Token::endOfLine;

      //std::cout << " " << token << std::flush;
    }



    std::ostream &operator<< ( std::ostream &out, const ProjectionBlock::Token &token )
    {
      typedef ProjectionBlock::Token Token;
      switch( token.type )
      {
      case Token::string:
        return out << "string [" << token.literal << "]";
      case Token::number:
        return out << "number [" << token.value << "]";
      case Token::defaultKeyword:
        return out << "default";
      case Token::functionKeyword:
        return out << "function";
      case Token::segmentKeyword:
        return out << "segment";
      case Token::sqrtKeyword:
        return out << "sqrt";
      case Token::sinKeyword:
        return out << "sin";
      case Token::cosKeyword:
        return out << "cos";
      case Token::piKeyword:
        return out << "pi";
      case Token::equals:
        return out << "'='";
      case Token::openingParen:
        return out << "'('";
      case Token::closingParen:
        return out << "')'";
      case Token::openingBracket:
        return out << "'['";
      case Token::closingBracket:
        return out << "']'";
      case Token::normDelim:
        return out << "'|'";
      case Token::additiveOperator:
        return out << "addop [" << token.symbol << "]";
      case Token::multiplicativeOperator:
        return out << "mulop [" << token.symbol << "]";
      case Token::powerOperator:
        return out << "powerop" << std::endl;
      case Token::endOfLine:
        return out << "eol";
      default:
        return out << "invalid [" << token.type << "]";
      }
    }


    void ProjectionBlock::ConstantExpression
      ::evaluate ( const Vector &argument, Vector &result ) const
    {
      result = value_;
    }


    void ProjectionBlock::VariableExpression
      ::evaluate ( const Vector &argument, Vector &result ) const
    {
      result = argument;
    }


    void ProjectionBlock::FunctionCallExpression
      ::evaluate ( const Vector &argument, Vector &result ) const
    {
      expression_->evaluate( argument, tmp_ );
      return function_->evaluate( tmp_, result );
    }


    void ProjectionBlock::BracketExpression
      ::evaluate ( const Vector &argument, Vector &result ) const
    {
      expression_->evaluate( argument, result );
      if( field_ >= result.size() )
        DUNE_THROW( MathError, "Index out of bounds (" <<  field_ << " not in [ 0, " << result.size() << " [)." );
      result[ 0 ] = result[ field_ ];
      result.resize( 1 );
    }


    void ProjectionBlock::MinusExpression
      ::evaluate ( const Vector &argument, Vector &result ) const
    {
      expression_->evaluate( argument, result );
      const size_t size = result.size();
      for( size_t i = 0; i < size; ++i )
        result[ i ] *= -1.0;
    }


    void ProjectionBlock::NormExpression
      ::evaluate ( const Vector &argument, Vector &result ) const
    {
      expression_->evaluate( argument, result );
      double normsqr = 0.0;
      const size_t size = result.size();
      for( size_t i = 0; i < size; ++i )
        normsqr += result[ i ] * result[ i ];
      result.resize( 1 );
      result[ 0 ] = sqrt( normsqr );
    }


    void ProjectionBlock::SqrtExpression
      ::evaluate ( const Vector &argument, Vector &result ) const
    {
      expression_->evaluate( argument, result );
      if( result.size() != 1 )
        DUNE_THROW( MathError, "Cannot calculate square root of a vector." );
      result[ 0 ] = sqrt( result[ 0 ] );
    }


    void ProjectionBlock::SinExpression
      ::evaluate ( const Vector &argument, Vector &result ) const
    {
      expression_->evaluate( argument, result );
      if( result.size() != 1 )
        DUNE_THROW( MathError, "Cannot calculate the sine of a vector." );
      result[ 0 ] = sin( result[ 0 ] );
    }


    void ProjectionBlock::CosExpression
      ::evaluate ( const Vector &argument, Vector &result ) const
    {
      expression_->evaluate( argument, result );
      if( result.size() != 1 )
        DUNE_THROW( MathError, "Cannot calculate the cosine of a vector." );
      result[ 0 ] = cos( result[ 0 ] );
    }


    void ProjectionBlock::PowerExpression
      ::evaluate ( const Vector &argument, Vector &result ) const
    {
      exprA_->evaluate( argument, result );
      exprB_->evaluate( argument, tmp_ );

      if( (result.size() == 1) && (tmp_.size() == 1) )
        result[ 0 ] = std::pow( result[ 0 ], tmp_[ 0 ] );
      else
        DUNE_THROW( MathError, "Cannot calculate powers of vectors." );
    }


    void ProjectionBlock::SumExpression
      ::evaluate ( const Vector &argument, Vector &result ) const
    {
      exprA_->evaluate( argument, result );
      exprB_->evaluate( argument, tmp_ );

      if( result.size() == tmp_.size() )
      {
        const size_t size = result.size();
        for( size_t i = 0; i < size; ++i )
          result[ i ] += tmp_[ i ];
      }
      else
        DUNE_THROW( MathError, "Cannot sum vectors of different size." );
    }


    void ProjectionBlock::DifferenceExpression
      ::evaluate ( const Vector &argument, Vector &result ) const
    {
      exprA_->evaluate( argument, result );
      exprB_->evaluate( argument, tmp_ );

      if( result.size() == tmp_.size() )
      {
        const size_t size = result.size();
        for( size_t i = 0; i < size; ++i )
          result[ i ] -= tmp_[ i ];
      }
      else
        DUNE_THROW( MathError, "Cannot sum vectors of different size." );
    }


    void ProjectionBlock::ProductExpression
      ::evaluate ( const Vector &argument, Vector &result ) const
    {
      exprA_->evaluate( argument, result );
      exprB_->evaluate( argument, tmp_ );

      if( result.size() == tmp_.size() )
      {
        double product = 0.0;
        const size_t size = result.size();
        for( size_t i = 0; i < size; ++i )
          product += result[ i ] * tmp_[ i ];
        result.resize( 1 );
        result[ 0 ] = product;
      }
      else if( tmp_.size() == 1 )
      {
        const size_t size = result.size();
        for( size_t i = 0; i < size; ++i )
          result[ i ] *= tmp_[ 0 ];
      }
      else if( result.size() == 1 )
      {
        std::swap( result, tmp_ );
        const size_t size = result.size();
        for( size_t i = 0; i < size; ++i )
          result[ i ] *= tmp_[ 0 ];
      }
      else
        DUNE_THROW( MathError, "Cannot multiply non-scalar vectors of different size." );
    }


    void ProjectionBlock::QuotientExpression
      ::evaluate ( const Vector &argument, Vector &result ) const
    {
      exprB_->evaluate( argument, result );
      if( result.size() != 1 )
        DUNE_THROW( MathError, "Cannot divide by a vector." );
      double factor = 1.0 / result[ 0 ];

      exprA_->evaluate( argument, result );
      const size_t size = result.size();
      for( size_t i = 0; i < size; ++i )
        result[ i ] *= factor;
    }

  }

}
